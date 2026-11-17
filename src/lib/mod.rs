//! Group BAM records into read groups by the value of a tag.

pub mod group;
pub mod raw;

use std::collections::hash_map::Entry;
use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs;
use std::io::BufWriter;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use bstr::{BString, ByteSlice};
use fgumi_raw_bam::{RawTagsView, aux_data_slice, update_string_tag};
use log::info;
use noodles::sam;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::program::tag::{
    COMMAND_LINE, NAME, PREVIOUS_PROGRAM_ID, VERSION,
};
use noodles_bgzf::io::writer::CompressionLevel;
use rayon::ThreadPool;
use rayon::prelude::*;
use tempfile::NamedTempFile;

use crate::group::{Plan, read_group_id, value_text};
use crate::raw::{Records, Sink, push_record, record_name};

/// What to group, where to read it from, and where to write it.
#[derive(Clone, Debug)]
pub struct Args {
    /// The BAM to read, or `-` for standard input.
    pub input: PathBuf,
    /// The BAM to write, or `-` for standard output.
    pub output: PathBuf,
    /// The tag whose values name the read groups.
    pub tag: Tag,
    /// What to do with records that lack the tag.
    pub untagged: Untagged,
    /// A value that records without the tag are grouped under instead.
    pub untagged_value: Option<String>,
    /// Threads for decompressing, grouping and compressing records.
    pub threads: NonZero<usize>,
    /// BGZF compression level of the output.
    pub compression_level: CompressionLevel,
    /// The command line recorded in the output's `@PG` line.
    pub command_line: String,
}

/// What to do with a record that lacks the tag.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, clap::ValueEnum)]
pub enum Untagged {
    /// Leave the record in its original read group.
    Keep,
    /// Omit the record from the output.
    Drop,
    /// Stop with an error at the first such record.
    #[default]
    Error,
}

/// Counts from one run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Summary {
    /// Records read.
    pub records: u64,
    /// Records written into a read group named by their tag value.
    pub grouped: u64,
    /// Records without the tag that kept their read group.
    pub kept: u64,
    /// Records without the tag that were dropped.
    pub dropped: u64,
    /// Read groups the tag values created.
    pub read_groups: usize,
}

/// The read group ID of each group key.
type GroupIds = HashMap<Vec<u8>, BString>;

/// The type byte under which `--untagged-value` keys its read group, apart from
/// every SAM tag type.
const UNTAGGED_VALUE: u8 = b'U';

/// Runs of records handed to each thread per batch, so a slow run does not
/// leave the other threads idle.
const RUNS_PER_THREAD: usize = 4;

/// A record's original read group and the raw value of the grouping tag.
struct Tags<'a> {
    original: Option<&'a [u8]>,
    value: Option<(u8, &'a [u8])>,
}

/// Find the record's read group and grouping tag in one walk over its tags.
fn tags(record: &[u8], tag: [u8; 2]) -> Result<Tags<'_>> {
    let mut found = Tags {
        original: None,
        value: None,
    };
    for entry in &RawTagsView::new(aux_data_slice(record)) {
        if entry.tag == *b"RG" && found.original.is_none() {
            if entry.type_byte != b'Z' {
                bail!("the RG tag of a record is not a string");
            }
            let id = entry.value_bytes;
            found.original = Some(id.strip_suffix(&[0]).unwrap_or(id));
        }
        if entry.tag == tag && found.value.is_none() {
            found.value = Some((entry.type_byte, entry.value_bytes));
        }
    }
    Ok(found)
}

/// Write into `key` what names the record's new read group: its original read
/// group and its tag value, or the untagged value, returning that value and its
/// type. Returns none when the record has neither, unless untagged records are
/// an error.
fn group_key<'a>(
    record: &[u8],
    tags: &Tags<'a>,
    args: &'a Args,
    key: &mut Vec<u8>,
) -> Result<Option<(u8, &'a [u8])>> {
    let (kind, value) = match (tags.value, &args.untagged_value) {
        (Some(value), _) => value,
        (None, Some(value)) => (UNTAGGED_VALUE, value.as_bytes()),
        (None, None) if args.untagged == Untagged::Error => bail!(
            "the record {} has no {} tag; pass --untagged keep or drop to allow it",
            record_name(record),
            tag_display(args.tag)
        ),
        (None, None) => return Ok(None),
    };
    key.clear();
    match tags.original {
        Some(original) => {
            key.push(1);
            key.extend_from_slice(original);
        }
        None => key.push(0),
    }
    key.push(0);
    key.push(kind);
    key.extend_from_slice(value);
    Ok(Some((kind, value)))
}

/// The text a group key's value names its read group by.
fn key_text(kind: u8, value: &[u8], args: &Args) -> Result<String> {
    if kind == UNTAGGED_VALUE {
        return Ok(value.to_str_lossy().into_owned());
    }
    value_text(kind, value).with_context(|| format!("in the {} tag", tag_display(args.tag)))
}

fn tag_display(tag: Tag) -> String {
    let bytes: [u8; 2] = tag.into();
    String::from_utf8_lossy(&bytes).into_owned()
}

/// Where the second pass reads the records from.
enum Source {
    /// The input file itself, read a second time.
    Input(PathBuf),
    /// A copy of a stream that cannot be read twice.
    Copy(NamedTempFile),
}

impl Source {
    fn path(&self) -> &Path {
        match self {
            Self::Input(path) => path,
            Self::Copy(file) => file.path(),
        }
    }
}

/// Whether `path` names a regular file, which can be read twice.
fn rereadable(path: &Path) -> bool {
    path != Path::new("-") && fs::metadata(path).is_ok_and(|metadata| metadata.is_file())
}

/// The read groups one run of records names: each group key new to the run,
/// with its original read group and text, and the original read groups of
/// records without one.
#[derive(Default)]
struct RunPlan {
    keys: Vec<(Vec<u8>, Option<BString>, String)>,
    untagged: BTreeSet<Option<BString>>,
}

fn plan_run(run: &[&[u8]], args: &Args, tag: [u8; 2]) -> Result<RunPlan> {
    let mut plan = RunPlan::default();
    let mut seen = HashSet::new();
    let mut key = Vec::new();
    for record in run {
        let tags = tags(record, tag)?;
        if let Some((kind, value)) = group_key(record, &tags, args, &mut key)? {
            if seen.insert(key.clone()) {
                let text = key_text(kind, value, args)?;
                plan.keys
                    .push((key.clone(), tags.original.map(BString::from), text));
            }
        } else {
            plan.untagged.insert(tags.original.map(BString::from));
        }
    }
    Ok(plan)
}

/// Learn every read group the tag implies, and the ID of each group key,
/// copying a stream as it goes so the records can be read again.
fn plan(args: &Args, pool: &ThreadPool) -> Result<(sam::Header, Plan, GroupIds, Source)> {
    let tag: [u8; 2] = args.tag.into();
    let (copy, source) = if rereadable(&args.input) {
        (None, Source::Input(args.input.clone()))
    } else {
        let file = tempfile::Builder::new()
            .prefix("bamby.")
            .suffix(".bam")
            .tempfile()
            .context("cannot create a temporary file for the input stream")?;
        let writer = BufWriter::new(file.reopen()?);
        (
            Some(Box::new(writer) as Box<dyn std::io::Write + Send>),
            Source::Copy(file),
        )
    };
    let (header, mut records) = Records::new(raw::open(&args.input)?, copy, pool)?;
    let mut plan = Plan::default();
    let mut ids = GroupIds::new();
    let parts = pool.current_num_threads() * RUNS_PER_THREAD;
    while let Some(batch) = records.next_batch(pool)? {
        let runs = batch.runs(parts);
        let planned: Vec<Result<RunPlan>> = pool.install(|| {
            runs.par_iter()
                .map(|run| plan_run(run, args, tag))
                .collect()
        });
        for run in planned {
            let run = run?;
            for (key, original, text) in run.keys {
                if let Entry::Vacant(slot) = ids.entry(key) {
                    let original = original.as_ref().map(|id| id.as_bstr());
                    slot.insert(read_group_id(original, &text));
                    plan.observe(original, Some(text));
                }
            }
            for original in run.untagged {
                plan.observe(original.as_ref().map(|id| id.as_bstr()), None);
            }
        }
    }
    Ok((header, plan, ids, source))
}

/// Add one `@PG` line for this run, following the last program in the header
/// as samtools does, with an ID made unique by a numeric suffix.
fn add_program(header: &mut sam::Header, command_line: &str) -> Result<()> {
    let programs = header.programs_mut().as_mut();
    let name = env!("CARGO_PKG_NAME");
    let mut id = BString::from(name);
    let mut suffix = 0;
    while programs.contains_key(&id) {
        suffix += 1;
        id = BString::from(format!("{name}.{suffix}"));
    }
    let mut builder = Map::<Program>::builder()
        .insert(NAME, name)
        .insert(VERSION, env!("CARGO_PKG_VERSION"))
        .insert(COMMAND_LINE, command_line);
    if let Some(previous) = programs.keys().last() {
        builder = builder.insert(PREVIOUS_PROGRAM_ID, previous.clone());
    }
    programs.insert(id, builder.build()?);
    Ok(())
}

/// Write every record into the read group of its tag value.
pub fn run(args: &Args) -> Result<Summary> {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads.get())
        .build()
        .context("cannot start the worker threads")?;
    let (input_header, plan, ids, source) = plan(args, &pool)?;
    let mut header = plan.header(&input_header, args.untagged == Untagged::Keep)?;
    add_program(&mut header, &args.command_line)?;

    let (_, mut records) = Records::new(raw::open(source.path())?, None, &pool)?;
    let level = u32::from(u8::from(args.compression_level));
    let mut sink = Sink::new(raw::create(&args.output)?, &header, level)?;

    let tag: [u8; 2] = args.tag.into();
    let parts = pool.current_num_threads() * RUNS_PER_THREAD;
    let mut summary = Summary {
        read_groups: plan.len(),
        ..Summary::default()
    };
    while let Some(batch) = records.next_batch(&pool)? {
        let runs = batch.runs(parts);
        let counts = sink.write_runs(&runs, &pool, |run, out| {
            let mut counts = Summary::default();
            let mut key = Vec::new();
            let mut edited = Vec::new();
            for record in run {
                counts.records += 1;
                let tags = tags(record, tag)?;
                if group_key(record, &tags, args, &mut key)?.is_some() {
                    let id = ids
                        .get(&key)
                        .context("a record names a read group the first pass did not see")?;
                    edited.clear();
                    edited.extend_from_slice(record);
                    update_string_tag(&mut edited, *b"RG", id);
                    push_record(out, &edited)?;
                    counts.grouped += 1;
                } else if args.untagged == Untagged::Drop {
                    counts.dropped += 1;
                } else {
                    push_record(out, record)?;
                    counts.kept += 1;
                }
            }
            Ok(counts)
        })?;
        for count in counts {
            summary.records += count.records;
            summary.grouped += count.grouped;
            summary.kept += count.kept;
            summary.dropped += count.dropped;
        }
    }
    sink.finish()?;

    info!(
        "Wrote {} of {} records into {} read groups; {} without the tag kept, {} dropped",
        summary.grouped, summary.records, summary.read_groups, summary.kept, summary.dropped
    );
    Ok(summary)
}
