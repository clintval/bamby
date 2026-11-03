//! Learn the read groups a tag implies and build the header that declares them.

use std::collections::{BTreeMap, BTreeSet};

use anyhow::{Result, bail};
use bstr::{BStr, BString, ByteSlice};
use noodles::sam;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::map::read_group::tag::SAMPLE;

/// The text of a tag value used to name a read group, or an error for a type
/// that has no single text form (floats and arrays).
pub fn value_text(value: &Value<'_>) -> Result<String> {
    let text = match value {
        Value::String(text) | Value::Hex(text) => text.to_str_lossy().into_owned(),
        Value::Character(c) => char::from(*c).to_string(),
        Value::Int8(n) => n.to_string(),
        Value::UInt8(n) => n.to_string(),
        Value::Int16(n) => n.to_string(),
        Value::UInt16(n) => n.to_string(),
        Value::Int32(n) => n.to_string(),
        Value::UInt32(n) => n.to_string(),
        Value::Float(_) => bail!("a float tag value cannot name a read group"),
        Value::Array(_) => bail!("an array tag value cannot name a read group"),
    };
    if text.is_empty() {
        bail!("an empty tag value cannot name a read group");
    }
    Ok(text)
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
    /// Record one record carrying `value` (or no value) under `original`.
    pub fn observe(&mut self, original: Option<&BStr>, value: Option<String>) {
        match value {
            Some(value) => {
                let id = read_group_id(original, &value);
                self.assignments.entry(id).or_insert_with(|| Assignment {
                    original: original.map(BString::from),
                    value,
                });
            }
            None => {
                self.untagged.insert(original.map(BString::from));
            }
        }
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
    #[case(Value::String(b"AAC-GTT".as_bstr()), "AAC-GTT")]
    #[case(Value::Hex(b"CAFE".as_bstr()), "CAFE")]
    #[case(Value::Character(b'x'), "x")]
    #[case(Value::Int8(-3), "-3")]
    #[case(Value::UInt16(7), "7")]
    #[case(Value::UInt32(4_000_000_000), "4000000000")]
    fn test_value_text_formats_scalar_values(#[case] value: Value<'_>, #[case] expected: &str) {
        assert_eq!(value_text(&value).unwrap(), expected);
    }

    #[test]
    fn test_value_text_rejects_floats_and_empty_strings() {
        assert!(value_text(&Value::Float(1.5)).is_err());
        assert!(value_text(&Value::String(b"".as_bstr())).is_err());
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
        plan.observe(Some(b"A".as_bstr()), Some("TTT".into()));
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()));
        plan.observe(Some(b"A".as_bstr()), Some("TTT".into()));

        let header = plan.header(&input, true).unwrap();

        assert_eq!(plan.len(), 2);
        assert_eq!(ids(&header), vec!["A.AAA", "A.TTT"]);
    }

    #[test]
    fn test_header_copies_the_original_fields_and_sets_the_sample() {
        let input = header_with_read_groups(&["A"]);
        let mut plan = Plan::default();
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()));

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
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()));
        plan.observe(Some(b"B".as_bstr()), None);

        assert_eq!(ids(&plan.header(&input, true).unwrap()), vec!["B", "A.AAA"]);
        assert_eq!(ids(&plan.header(&input, false).unwrap()), vec!["A.AAA"]);
    }

    #[test]
    fn test_header_names_values_of_records_without_a_read_group_by_value_alone() {
        let input = sam::Header::default();
        let mut plan = Plan::default();
        plan.observe(None, Some("AAA".into()));

        let header = plan.header(&input, true).unwrap();
        let map = header.read_groups().get(b"AAA".as_bstr()).unwrap();

        assert_eq!(map.other_fields().get(&SAMPLE).unwrap(), "AAA");
        assert_eq!(map.other_fields().get(&LIBRARY), None);
    }

    #[test]
    fn test_header_rejects_a_new_read_group_that_collides_with_a_kept_one() {
        let input = header_with_read_groups(&["A.AAA", "A"]);
        let mut plan = Plan::default();
        plan.observe(Some(b"A".as_bstr()), Some("AAA".into()));
        plan.observe(Some(b"A.AAA".as_bstr()), None);

        let error = plan.header(&input, true).unwrap_err().to_string();

        assert!(error.contains("collides"), "got: {error}");
    }
}
