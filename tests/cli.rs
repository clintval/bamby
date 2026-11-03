//! The real binary run on small BAMs built in a temporary directory.

use std::fs::File;
use std::io::Cursor;
use std::path::{Path, PathBuf};

use assert_cmd::Command;
use bstr::{BStr, ByteSlice};
use noodles::bam;
use noodles::sam;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Data, QualityScores, RecordBuf, Sequence};
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::program::tag::{NAME, PREVIOUS_PROGRAM_ID};
use noodles::sam::header::record::value::map::read_group::tag::{LIBRARY, PLATFORM, SAMPLE};
use noodles::sam::header::record::value::map::{Program, ReadGroup};
use predicates::str::contains;
use pretty_assertions::assert_eq;
use tempfile::TempDir;

const CB: Tag = Tag::new(b'C', b'B');

/// An unmapped record named `name` with the given tags.
fn record(name: &str, tags: Vec<(Tag, Value)>) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::UNMAPPED)
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_quality_scores(QualityScores::from(vec![30, 30, 30, 30]))
        .set_data(tags.into_iter().collect::<Data>())
        .build()
}

/// A record in read group `A` with cell barcode `cb`, if any.
fn cell(name: &str, cb: Option<&str>) -> RecordBuf {
    let mut tags = vec![(Tag::READ_GROUP, Value::from("A"))];
    if let Some(cb) = cb {
        tags.push((CB, Value::from(cb)));
    }
    record(name, tags)
}

/// A header with read group `A` from library `lib`.
fn header() -> sam::Header {
    let map = Map::<ReadGroup>::builder()
        .insert(SAMPLE, "pool")
        .insert(LIBRARY, "lib")
        .insert(PLATFORM, "ILLUMINA")
        .build()
        .unwrap();
    sam::Header::builder().add_read_group("A", map).build()
}

fn write_bam(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    let mut writer = bam::io::Writer::new(File::create(path).unwrap());
    writer.write_header(header).unwrap();
    for record in records {
        writer.write_alignment_record(header, record).unwrap();
    }
    writer.try_finish().unwrap();
}

fn read_bam<R: std::io::Read>(inner: R) -> (sam::Header, Vec<RecordBuf>) {
    let mut reader = bam::io::Reader::new(inner);
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .collect::<Result<Vec<_>, _>>()
        .unwrap();
    (header, records)
}

/// A temporary directory holding `in.bam` with the four standard records.
fn fixture() -> (TempDir, PathBuf) {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    let records = [
        cell("r1", Some("TTT")),
        cell("r2", Some("AAA")),
        cell("r3", Some("TTT")),
        cell("r4", None),
    ];
    write_bam(&input, &header(), &records);
    (dir, input)
}

fn bamby() -> Command {
    let mut cmd = Command::cargo_bin("bamby").unwrap();
    cmd.env("NO_COLOR", "1");
    cmd
}

fn read_group_ids(header: &sam::Header) -> Vec<String> {
    header
        .read_groups()
        .keys()
        .map(ToString::to_string)
        .collect()
}

fn names_and_read_groups(records: &[RecordBuf]) -> Vec<(String, String)> {
    records
        .iter()
        .map(|record| {
            let name = record.name().unwrap().to_string();
            let rg = match record.data().get(&Tag::READ_GROUP) {
                Some(Value::String(id)) => id.to_string(),
                other => panic!("unexpected RG: {other:?}"),
            };
            (name, rg)
        })
        .collect()
}

fn pairs(expected: &[(&str, &str)]) -> Vec<(String, String)> {
    expected
        .iter()
        .map(|(a, b)| ((*a).to_string(), (*b).to_string()))
        .collect()
}

#[test]
fn test_every_value_gets_its_own_read_group() {
    let (dir, input) = fixture();
    let output = dir.path().join("out.bam");

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["-t", "CB"])
        .assert()
        .success();

    let (header, records) = read_bam(File::open(&output).unwrap());
    assert_eq!(read_group_ids(&header), vec!["A", "A.AAA", "A.TTT"]);
    assert_eq!(
        names_and_read_groups(&records),
        pairs(&[
            ("r1", "A.TTT"),
            ("r2", "A.AAA"),
            ("r3", "A.TTT"),
            ("r4", "A")
        ])
    );
}

#[test]
fn test_new_read_groups_copy_the_original_fields_and_take_the_value_as_sample() {
    let (dir, input) = fixture();
    let output = dir.path().join("out.bam");

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["-t", "CB"])
        .assert()
        .success();

    let (header, _) = read_bam(File::open(&output).unwrap());
    let map = header.read_groups().get(b"A.TTT".as_bstr()).unwrap();
    assert_eq!(map.other_fields().get(&SAMPLE).unwrap(), "TTT");
    assert_eq!(map.other_fields().get(&LIBRARY).unwrap(), "lib");
    assert_eq!(map.other_fields().get(&PLATFORM).unwrap(), "ILLUMINA");
}

#[test]
fn test_drop_untagged_removes_the_records_and_their_unused_read_group() {
    let (dir, input) = fixture();
    let output = dir.path().join("out.bam");

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["-t", "CB", "--drop-untagged"])
        .assert()
        .success();

    let (header, records) = read_bam(File::open(&output).unwrap());
    assert_eq!(read_group_ids(&header), vec!["A.AAA", "A.TTT"]);
    assert_eq!(
        names_and_read_groups(&records),
        pairs(&[("r1", "A.TTT"), ("r2", "A.AAA"), ("r3", "A.TTT")])
    );
}

#[test]
fn test_records_without_a_read_group_are_named_by_value_alone() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let xy = Tag::new(b'X', b'Y');
    write_bam(
        &input,
        &sam::Header::default(),
        &[record("r1", vec![(xy, Value::from(7_u8))])],
    );

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["--tag", "XY"])
        .assert()
        .success();

    let (header, records) = read_bam(File::open(&output).unwrap());
    assert_eq!(read_group_ids(&header), vec!["7"]);
    assert_eq!(names_and_read_groups(&records), pairs(&[("r1", "7")]));
}

#[test]
fn test_output_to_standard_output() {
    let (_dir, input) = fixture();

    let assert = bamby()
        .arg(&input)
        .args(["-o", "-", "-t", "CB", "-l", "0"])
        .assert()
        .success();

    let (header, records) = read_bam(Cursor::new(assert.get_output().stdout.clone()));
    assert_eq!(read_group_ids(&header), vec!["A", "A.AAA", "A.TTT"]);
    assert_eq!(records.len(), 4);
}

#[test]
fn test_threads_write_the_same_records() {
    let (dir, input) = fixture();
    let single = dir.path().join("single.bam");
    let multi = dir.path().join("multi.bam");

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&single)
        .args(["-t", "CB"])
        .assert()
        .success();
    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&multi)
        .args(["-t", "CB", "-@", "4"])
        .assert()
        .success();

    let (_, single) = read_bam(File::open(&single).unwrap());
    let (_, multi) = read_bam(File::open(&multi).unwrap());
    assert_eq!(single, multi);
}

#[test]
fn test_the_run_is_recorded_as_one_program_following_the_last() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    let program = |name: &str| {
        Map::<Program>::builder()
            .insert(NAME, name)
            .build()
            .unwrap()
    };
    let header = sam::Header::builder()
        .add_program("aligner", program("aligner"))
        .add_program("other", program("other"))
        .add_program("bamby", program("bamby"))
        .build();
    write_bam(
        &input,
        &header,
        &[record("r1", vec![(CB, Value::from("AAA"))])],
    );

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["-t", "CB"])
        .assert()
        .success();

    let (header, _) = read_bam(File::open(&output).unwrap());
    let programs = header.programs().as_ref();
    let ids: Vec<String> = programs.keys().map(ToString::to_string).collect();
    assert_eq!(ids, vec!["aligner", "other", "bamby", "bamby.1"]);
    let ours = programs.get(b"bamby.1".as_bstr()).unwrap();
    assert_eq!(
        ours.other_fields().get(&NAME).map(|n| n.as_bstr()),
        Some(BStr::new("bamby"))
    );
    assert_eq!(
        ours.other_fields()
            .get(&PREVIOUS_PROGRAM_ID)
            .map(|n| n.as_bstr()),
        Some(BStr::new("bamby"))
    );
}

#[test]
fn test_standard_input_is_copied_and_written_to_standard_output() {
    let (dir, input) = fixture();
    let output = dir.path().join("out.bam");

    bamby()
        .arg(&input)
        .args(["-t", "CB", "-o"])
        .arg(&output)
        .assert()
        .success();
    let assert = bamby()
        .args(["-t", "CB"])
        .pipe_stdin(&input)
        .unwrap()
        .assert()
        .success();

    let (_, from_file) = read_bam(File::open(&output).unwrap());
    let (header, from_stdin) = read_bam(Cursor::new(assert.get_output().stdout.clone()));
    assert_eq!(read_group_ids(&header), vec!["A", "A.AAA", "A.TTT"]);
    assert_eq!(from_stdin, from_file);
}

#[test]
fn test_a_float_tag_is_an_error() {
    let dir = TempDir::new().unwrap();
    let input = dir.path().join("in.bam");
    let output = dir.path().join("out.bam");
    write_bam(
        &input,
        &sam::Header::default(),
        &[record("r1", vec![(CB, Value::from(1.5_f32))])],
    );

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .args(["-t", "CB"])
        .assert()
        .failure()
        .stderr(contains("float tag value"));
}

#[test]
fn test_a_missing_input_is_an_error() {
    let dir = TempDir::new().unwrap();

    bamby()
        .arg(dir.path().join("missing.bam"))
        .arg("-o")
        .arg(dir.path().join("out.bam"))
        .args(["-t", "CB"])
        .assert()
        .failure()
        .stderr(contains("cannot open"));
}

#[test]
fn test_the_tag_is_required() {
    let (dir, input) = fixture();

    bamby()
        .arg(&input)
        .arg("-o")
        .arg(dir.path().join("out.bam"))
        .assert()
        .failure()
        .stderr(contains("--tag <TAG>"));
}
