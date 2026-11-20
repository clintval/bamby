//! Learn the read groups a tag implies and build the header that declares them.

use std::collections::{BTreeMap, BTreeSet};

use anyhow::{Result, bail};
use bstr::{BStr, BString, ByteSlice};
use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag::SAMPLE;

/// The text of a raw tag value, of type `kind`, used to name a read group, or
/// an error for a type that has no single text form (floats and arrays).
pub fn value_text(kind: u8, bytes: &[u8]) -> Result<String> {
    let text = match kind {
        b'Z' | b'H' => bytes
            .strip_suffix(&[0])
            .unwrap_or(bytes)
            .to_str_lossy()
            .into_owned(),
        b'A' => char::from(fixed::<1>(bytes)?[0]).to_string(),
        b'c' => i8::from_le_bytes(fixed(bytes)?).to_string(),
        b'C' => u8::from_le_bytes(fixed(bytes)?).to_string(),
        b's' => i16::from_le_bytes(fixed(bytes)?).to_string(),
        b'S' => u16::from_le_bytes(fixed(bytes)?).to_string(),
        b'i' => i32::from_le_bytes(fixed(bytes)?).to_string(),
        b'I' => u32::from_le_bytes(fixed(bytes)?).to_string(),
        b'f' => bail!("a float tag value cannot name a read group"),
        b'B' => bail!("an array tag value cannot name a read group"),
        _ => bail!(
            "a tag of type {} cannot name a read group",
            char::from(kind)
        ),
    };
    if text.is_empty() {
        bail!("an empty tag value cannot name a read group");
    }
    Ok(text)
}

fn fixed<const N: usize>(bytes: &[u8]) -> Result<[u8; N]> {
    match bytes.get(..N) {
        Some(head) => Ok(head.try_into()?),
        None => bail!("a tag value is shorter than its type"),
    }
}

/// The ID of the read group for `value`: the record's original read group and
/// the value joined by a period, or the value alone when it has none.
pub fn read_group_id(original: Option<&BStr>, value: &str) -> BString {
    match original {
        Some(original) => {
            let mut id = BString::from(original.to_vec());
            id.push(b'.');
            id.extend_from_slice(value.as_bytes());
            id
        }
        None => BString::from(value),
    }
}

/// One new read group: the value it holds and the read group it came from.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Assignment {
    original: Option<BString>,
    value: String,
}

/// The read groups seen in a first pass over the records.
#[derive(Debug, Default)]
pub struct Plan {
    assignments: BTreeMap<BString, Assignment>,
    untagged: BTreeSet<Option<BString>>,
}

impl Plan {
    /// Record one record carrying `value` (or no value) under `original`, or an
    /// error when another value already has the same read group ID.
    pub fn observe(&mut self, original: Option<&BStr>, value: Option<String>) -> Result<()> {
        let Some(value) = value else {
            self.untagged.insert(original.map(BString::from));
            return Ok(());
        };
        let id = read_group_id(original, &value);
        let assignment = Assignment {
            original: original.map(BString::from),
            value,
        };
        let existing = self
            .assignments
            .entry(id.clone())
            .or_insert_with(|| assignment.clone());
        if *existing != assignment {
            bail!(
                "the tag values {} and {} would share the read group {id}",
                existing.value,
                assignment.value
            );
        }
        Ok(())
    }

    /// The number of read groups the tag values create.
    pub fn len(&self) -> usize {
        self.assignments.len()
    }

    /// Whether no record carried the tag.
    pub fn is_empty(&self) -> bool {
        self.assignments.is_empty()
    }

    /// The output header: `input` with one read group per value, and with only
    /// those original read groups that untagged records still use when they are kept.
    pub fn header(&self, input: &sam::Header, keep_untagged: bool) -> Result<sam::Header> {
        let mut header = input.clone();
        let originals = input.read_groups();

        let retained: BTreeSet<BString> = if keep_untagged {
            self.untagged.iter().flatten().cloned().collect()
        } else {
            BTreeSet::new()
        };

        let read_groups = header.read_groups_mut();
        read_groups.retain(|id, _| retained.contains(id));

        for (id, assignment) in &self.assignments {
            if read_groups.contains_key(id) {
                bail!(
                    "the read group {id} for tag value {} collides with an existing read group",
                    assignment.value
                );
            }
            let mut map = assignment
                .original
                .as_ref()
                .and_then(|original| originals.get(original))
                .cloned()
                .unwrap_or_else(Map::<ReadGroup>::default);
            map.other_fields_mut()
                .insert(SAMPLE, BString::from(assignment.value.as_str()));
            read_groups.insert(id.clone(), map);
        }

        Ok(header)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use noodles::sam::header::record::value::map::read_group::tag::{LIBRARY, PLATFORM};
    use pretty_assertions::assert_eq;
    use rstest::rstest;

    fn header_with_read_groups(ids: &[&str]) -> sam::Header {
        let mut builder = sam::Header::builder();
        for id in ids {
            let map = Map::<ReadGroup>::builder()
                .insert(SAMPLE, format!("{id}-sample"))
                .insert(LIBRARY, "lib")
                .insert(PLATFORM, "ILLUMINA")
                .build()
                .unwrap();
            builder = builder.add_read_group(*id, map);
        }
        builder.build()
    }

    fn ids(header: &sam::Header) -> Vec<String> {
        header
            .read_groups()
            .keys()
            .map(ToString::to_string)
            .collect()
    }

    #[rstest]
    #[case(b'Z', b"AAC-GTT\0".to_vec(), "AAC-GTT")]
    #[case(b'H', b"CAFE\0".to_vec(), "CAFE")]
    #[case(b'A', b"x".to_vec(), "x")]
    #[case(b'c', (-3i8).to_le_bytes().to_vec(), "-3")]
    #[case(b'S', 7u16.to_le_bytes().to_vec(), "7")]
    #[case(b'I', 4_000_000_000u32.to_le_bytes().to_vec(), "4000000000")]
    fn test_value_text_formats_scalar_values(
        #[case] kind: u8,
        #[case] bytes: Vec<u8>,
        #[case] expected: &str,
    ) {
        assert_eq!(value_text(kind, &bytes).unwrap(), expected);
    }

    #[test]
    fn test_value_text_rejects_floats_arrays_and_empty_strings() {
        assert!(value_text(b'f', &1.5f32.to_le_bytes()).is_err());
        assert!(value_text(b'B', b"C\x01\x00\x00\x00\x07").is_err());
        assert!(value_text(b'Z', b"\0").is_err());
    }

    #[test]
    fn test_read_group_id_joins_the_original_and_the_value() {
        assert_eq!(read_group_id(Some(b"A".as_bstr()), "AAC"), "A.AAC");
        assert_eq!(read_group_id(None, "AAC"), "AAC");
    }

    #[test]
    fn test_header_declares_one_read_group_per_value_in_sorted_order() {
        let input = header_with_read_groups(&["A"]);
        let mut plan = Plan::default();
        plan.observe(Some(b"A".as_bstr()), Some("TTT".into()))
            .unwrap();
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()))
            .unwrap();
        plan.observe(Some(b"A".as_bstr()), Some("TTT".into()))
            .unwrap();

        let header = plan.header(&input, true).unwrap();

        assert_eq!(plan.len(), 2);
        assert_eq!(ids(&header), vec!["A.AAA", "A.TTT"]);
    }

    #[test]
    fn test_header_copies_the_original_fields_and_sets_the_sample() {
        let input = header_with_read_groups(&["A"]);
        let mut plan = Plan::default();
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()))
            .unwrap();

        let header = plan.header(&input, true).unwrap();
        let map = header.read_groups().get(b"A.AAA".as_bstr()).unwrap();

        assert_eq!(map.other_fields().get(&SAMPLE).unwrap(), "AAA");
        assert_eq!(map.other_fields().get(&LIBRARY).unwrap(), "lib");
        assert_eq!(map.other_fields().get(&PLATFORM).unwrap(), "ILLUMINA");
    }

    #[test]
    fn test_header_keeps_an_original_read_group_only_while_untagged_records_use_it() {
        let input = header_with_read_groups(&["A", "B"]);
        let mut plan = Plan::default();
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()))
            .unwrap();
        plan.observe(Some(b"B".as_bstr()), None).unwrap();

        assert_eq!(ids(&plan.header(&input, true).unwrap()), vec!["B", "A.AAA"]);
        assert_eq!(ids(&plan.header(&input, false).unwrap()), vec!["A.AAA"]);
    }

    #[test]
    fn test_header_names_values_of_records_without_a_read_group_by_value_alone() {
        let input = sam::Header::default();
        let mut plan = Plan::default();
        plan.observe(None, Some("AAA".into())).unwrap();

        let header = plan.header(&input, true).unwrap();
        let map = header.read_groups().get(b"AAA".as_bstr()).unwrap();

        assert_eq!(map.other_fields().get(&SAMPLE).unwrap(), "AAA");
        assert_eq!(map.other_fields().get(&LIBRARY), None);
    }

    #[test]
    fn test_header_rejects_a_new_read_group_that_collides_with_a_kept_one() {
        let input = header_with_read_groups(&["A.AAA", "A"]);
        let mut plan = Plan::default();
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()))
            .unwrap();
        plan.observe(Some(b"A.AAA".as_bstr()), None).unwrap();

        let error = plan.header(&input, true).unwrap_err().to_string();

        assert!(error.contains("collides"), "got: {error}");
    }

    #[test]
    fn test_observe_rejects_two_values_that_would_share_a_read_group() {
        let mut plan = Plan::default();
        plan.observe(Some(b"A".as_bstr()), Some("B.C".into()))
            .unwrap();

        let error = plan
            .observe(Some(b"A.B".as_bstr()), Some("C".into()))
            .unwrap_err();

        assert!(
            error.to_string().contains("share the read group A.B.C"),
            "got: {error}"
        );
    }
}
