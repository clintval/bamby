//! Group BAM records into read groups by the value of a tag.

pub mod group;
pub mod io;

use std::fs;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use bstr::{BStr, BString};
use log::info;
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::data::field::{Tag, Value};
use noodles::sam::alignment::record_buf::RecordBuf;
use noodles::sam::alignment::record_buf::data::field::Value as BufValue;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::program::tag::{
    COMMAND_LINE, NAME, PREVIOUS_PROGRAM_ID, VERSION,
};
use noodles_bgzf::io::writer::CompressionLevel;
use tempfile::NamedTempFile;

use crate::group::{Plan, read_group_id, value_text};

/// What to group, where to read it from, and where to write it.
#[derive(Clone, Debug)]
pub struct Args {
    /// The BAM to read, or `-` for standard input.
    pub input: PathBuf,
    /// The BAM to write, or `-` for standard output.
    pub output: PathBuf,
    /// The tag whose values name the read groups.
    pub tag: Tag,
    /// Drop records without the tag instead of keeping their read group.
    pub drop_untagged: bool,
    /// BGZF worker threads for reading and for writing.
    pub threads: NonZero<usize>,
    /// BGZF compression level of the output.
    pub compression_level: CompressionLevel,
    /// The command line recorded in the output's `@PG` line.
    pub command_line: String,
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

/// The record's read group, borrowed from its data.
fn original_read_group(record: &bam::Record) -> Result<Option<&BStr>> {
    match record.data().get(&Tag::READ_GROUP).transpose()? {
        Some(Value::String(id)) => Ok(Some(id)),
        Some(_) => bail!("the RG tag of a record is not a string"),
        None => Ok(None),
    }
}

/// The text of the record's `tag`, if it has one.
fn tag_text(record: &bam::Record, tag: Tag) -> Result<Option<String>> {
    match record.data().get(&tag).transpose()? {
        Some(value) => value_text(&value)
            .map(Some)
            .with_context(|| format!("in the {} tag", tag_display(tag))),
        None => Ok(None),
    }
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

/// Learn every read group the tag implies, copying a stream as it goes so the
/// records can be read again.
fn plan(args: &Args) -> Result<(sam::Header, Plan, Source)> {
    let mut reader = io::open(&args.input, args.threads)?;
    let header = reader.read_header()?;
    let mut copy = if rereadable(&args.input) {
        None
    } else {
        let file = tempfile::Builder::new()
            .prefix("bamby.")
            .suffix(".bam")
            .tempfile()
            .context("cannot create a temporary file for the input stream")?;
        let sink = io::Sink::create(file.path(), args.threads, CompressionLevel::FAST)?;
        let mut writer = io::Writer::from(sink);
        writer.write_header(&header)?;
        Some((file, writer))
    };
    let mut plan = Plan::default();
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        plan.observe(original_read_group(&record)?, tag_text(&record, args.tag)?);
        if let Some((_, writer)) = copy.as_mut() {
            writer.write_record(&header, &record)?;
        }
    }
    let source = match copy {
        Some((file, writer)) => {
            writer.into_inner().finish()?;
            Source::Copy(file)
        }
        None => Source::Input(args.input.clone()),
    };
    Ok((header, plan, source))
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
    let (input_header, plan, source) = plan(args)?;
    let mut header = plan.header(&input_header, !args.drop_untagged)?;
    add_program(&mut header, &args.command_line)?;

    let mut reader = io::open(source.path(), args.threads)?;
    reader.read_header()?;
    let sink = io::Sink::create(&args.output, args.threads, args.compression_level)?;
    let mut writer = io::Writer::from(sink);
    writer.write_header(&header)?;

    let mut summary = Summary {
        read_groups: plan.len(),
        ..Summary::default()
    };
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        summary.records += 1;
        match tag_text(&record, args.tag)? {
            Some(value) => {
                let id: BString = read_group_id(original_read_group(&record)?, &value);
                let mut buf = RecordBuf::try_from_alignment_record(&input_header, &record)?;
                buf.data_mut().insert(Tag::READ_GROUP, BufValue::String(id));
                writer.write_alignment_record(&header, &buf)?;
                summary.grouped += 1;
            }
            None if args.drop_untagged => summary.dropped += 1,
            None => {
                writer.write_alignment_record(&header, &record)?;
                summary.kept += 1;
            }
        }
    }
    writer.into_inner().finish()?;

    info!(
        "Wrote {} of {} records into {} read groups; {} without the tag kept, {} dropped",
        summary.grouped, summary.records, summary.read_groups, summary.kept, summary.dropped
    );
    Ok(summary)
}
