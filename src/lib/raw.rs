//! BAM records as raw bytes: BGZF blocks read and decompressed in parallel, and
//! output compressed in parallel, with no record decoded into fields.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use fgumi_bgzf::{
    BGZF_EOF, Decompressor, InlineBgzfCompressor, RawBgzfBlock, decompress_block_into,
    read_raw_blocks,
};
use noodles::bam;
use noodles::sam;
use rayon::ThreadPool;
use rayon::prelude::*;

/// BGZF blocks read at a time: up to 16 MiB of records once decompressed.
const BLOCKS_PER_BATCH: usize = 256;

/// The BAM magic number that opens every BAM stream.
const BAM_MAGIC: &[u8; 4] = b"BAM\x01";

/// Open `path`, or standard input when it is `-`.
pub fn open(path: &Path) -> Result<Box<dyn Read + Send>> {
    if path == Path::new("-") {
        return Ok(Box::new(BufReader::new(io::stdin())));
    }
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    Ok(Box::new(BufReader::new(file)))
}

/// Create `path`, or standard output when it is `-`.
pub fn create(path: &Path) -> Result<Box<dyn Write + Send>> {
    if path == Path::new("-") {
        return Ok(Box::new(BufWriter::new(io::stdout())));
    }
    let file = File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
    Ok(Box::new(BufWriter::new(file)))
}

/// A BAM stream read as its header and then batches of whole records.
pub struct Records {
    reader: Box<dyn Read + Send>,
    copy: Option<Box<dyn Write + Send>>,
    pending: Vec<u8>,
    start: usize,
    exhausted: bool,
}

impl Records {
    /// Read the header of the BAM in `reader`, copying every compressed block
    /// to `copy`, when given, so the stream can be read again.
    pub fn new(
        reader: Box<dyn Read + Send>,
        copy: Option<Box<dyn Write + Send>>,
        pool: &ThreadPool,
    ) -> Result<(sam::Header, Self)> {
        let mut records = Self {
            reader,
            copy,
            pending: Vec::new(),
            start: 0,
            exhausted: false,
        };
        let length = loop {
            if let Some(length) = header_length(&records.pending)? {
                break length;
            }
            if !records.fill(pool)? {
                bail!("the input ends before its BAM header does");
            }
        };
        let header = bam::io::Reader::from(&records.pending[..length])
            .read_header()
            .context("cannot read the BAM header")?;
        records.start = length;
        Ok((header, records))
    }

    /// The next batch of whole records, each without its leading length, or
    /// none at the end of the stream.
    pub fn next_batch(&mut self, pool: &ThreadPool) -> Result<Option<Batch<'_>>> {
        self.pending.drain(..self.start);
        self.start = 0;
        let ends = loop {
            let ends = record_ends(&self.pending)?;
            if !ends.is_empty() {
                break ends;
            }
            if !self.fill(pool)? {
                if self.pending.is_empty() {
                    return Ok(None);
                }
                bail!("the input ends inside a BAM record");
            }
        };
        self.start = ends.last().copied().unwrap_or_default();
        Ok(Some(Batch {
            bytes: &self.pending[..self.start],
            ends,
        }))
    }

    /// Decompress the next blocks onto the pending bytes, returning false at
    /// the end of the stream.
    fn fill(&mut self, pool: &ThreadPool) -> Result<bool> {
        if self.exhausted {
            return Ok(false);
        }
        let blocks = read_raw_blocks(&mut self.reader, BLOCKS_PER_BATCH)
            .context("cannot read a BGZF block")?;
        if blocks.is_empty() {
            self.exhausted = true;
            if let Some(mut copy) = self.copy.take() {
                copy.write_all(&BGZF_EOF)?;
                copy.flush()?;
            }
            return Ok(false);
        }
        if let Some(copy) = self.copy.as_mut() {
            for block in &blocks {
                copy.write_all(&block.data)?;
            }
        }
        for bytes in inflate(&blocks, pool)? {
            self.pending.extend_from_slice(&bytes);
        }
        Ok(true)
    }
}

/// Whole records from one stretch of a BAM stream.
pub struct Batch<'a> {
    bytes: &'a [u8],
    ends: Vec<usize>,
}

impl<'a> Batch<'a> {
    /// The records, each without its leading length, split into about `parts`
    /// runs of neighbouring records, in order.
    pub fn runs(&self, parts: usize) -> Vec<Vec<&'a [u8]>> {
        let records: Vec<&'a [u8]> = self
            .ends
            .iter()
            .scan(0, |start, &end| {
                let record = &self.bytes[*start + 4..end];
                *start = end;
                Some(record)
            })
            .collect();
        let size = records.len().div_ceil(parts.max(1)).max(1);
        records.chunks(size).map(<[&[u8]]>::to_vec).collect()
    }
}

/// Decompress `blocks` in parallel, in order.
fn inflate(blocks: &[RawBgzfBlock], pool: &ThreadPool) -> Result<Vec<Vec<u8>>> {
    pool.install(|| {
        blocks
            .par_iter()
            .map_init(Decompressor::new, |decompressor, block| {
                let mut bytes = Vec::with_capacity(block.uncompressed_size());
                decompress_block_into(block, decompressor, &mut bytes).map(|()| bytes)
            })
            .collect::<io::Result<Vec<_>>>()
            .context("cannot decompress a BGZF block")
    })
}

/// The length of the BAM header at the start of `bytes`, or none until all of
/// it is there.
fn header_length(bytes: &[u8]) -> Result<Option<usize>> {
    if bytes.len() < 4 {
        return Ok(None);
    }
    if &bytes[..4] != BAM_MAGIC {
        bail!("the input is not a BAM file");
    }
    let mut at = 4;
    let Some(text) = read_u32(bytes, at) else {
        return Ok(None);
    };
    at += 4 + text as usize;
    let Some(references) = read_u32(bytes, at) else {
        return Ok(None);
    };
    at += 4;
    for _ in 0..references {
        let Some(name) = read_u32(bytes, at) else {
            return Ok(None);
        };
        at += 4 + name as usize + 4;
    }
    Ok((at <= bytes.len()).then_some(at))
}

/// The end of each whole record in `bytes`, after its leading length.
fn record_ends(bytes: &[u8]) -> Result<Vec<usize>> {
    let mut ends = Vec::new();
    let mut at = 0;
    while let Some(length) = read_u32(bytes, at) {
        if length < 32 {
            bail!("a BAM record of {length} bytes is too short to hold its fields");
        }
        let end = at + 4 + length as usize;
        if end > bytes.len() {
            break;
        }
        ends.push(end);
        at = end;
    }
    Ok(ends)
}

fn read_u32(bytes: &[u8], at: usize) -> Option<u32> {
    let field = bytes.get(at..at + 4)?;
    Some(u32::from_le_bytes(field.try_into().ok()?))
}

/// A BAM output: the header, then runs of records compressed in parallel.
pub struct Sink {
    writer: Box<dyn Write + Send>,
    level: u32,
}

impl Sink {
    /// Write `header` to `writer` as the start of a BAM at compression `level`.
    pub fn new(
        mut writer: Box<dyn Write + Send>,
        header: &sam::Header,
        level: u32,
    ) -> Result<Self> {
        let mut bytes = bam::io::Writer::from(Vec::new());
        bytes.write_header(header)?;
        let mut compressor = InlineBgzfCompressor::new(level);
        compressor.write_all(&bytes.into_inner())?;
        compressor.flush()?;
        compressor.write_blocks_to(&mut writer)?;
        Ok(Self { writer, level })
    }

    /// Fill one output run per input run with `fill` and compress each in
    /// parallel, then write them in order, returning what `fill` returned.
    pub fn write_runs<'a, T, F>(
        &mut self,
        runs: &[Vec<&'a [u8]>],
        pool: &ThreadPool,
        fill: F,
    ) -> Result<Vec<T>>
    where
        T: Send,
        F: Fn(&[&'a [u8]], &mut Vec<u8>) -> Result<T> + Sync,
    {
        let level = self.level;
        let compressed: Vec<(Vec<u8>, T)> = pool.install(|| {
            runs.par_iter()
                .map(|run| {
                    let mut out = Vec::new();
                    let value = fill(run, &mut out)?;
                    let mut compressor = InlineBgzfCompressor::new(level);
                    compressor.write_all(&out)?;
                    compressor.flush()?;
                    let mut blocks = Vec::new();
                    compressor.write_blocks_to(&mut blocks)?;
                    Ok((blocks, value))
                })
                .collect::<Result<Vec<_>>>()
        })?;
        let mut values = Vec::with_capacity(compressed.len());
        for (blocks, value) in compressed {
            self.writer.write_all(&blocks)?;
            values.push(value);
        }
        Ok(values)
    }

    /// End the BAM with the BGZF end-of-file marker.
    pub fn finish(mut self) -> Result<()> {
        self.writer.write_all(&BGZF_EOF)?;
        self.writer.flush()?;
        Ok(())
    }
}

/// Write `record` to `out` with its leading length.
pub fn push_record(out: &mut Vec<u8>, record: &[u8]) -> Result<()> {
    let length = u32::try_from(record.len()).context("a BAM record is longer than 4 GiB")?;
    out.extend_from_slice(&length.to_le_bytes());
    out.extend_from_slice(record);
    Ok(())
}

/// The name of `record`, as SAM writes it.
pub fn record_name(record: &[u8]) -> String {
    let length = usize::from(record[8]);
    let name = record
        .get(32..32 + length.saturating_sub(1))
        .unwrap_or_default();
    String::from_utf8_lossy(name).into_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    use pretty_assertions::assert_eq;

    fn record(length: u32) -> Vec<u8> {
        let mut bytes = length.to_le_bytes().to_vec();
        bytes.extend(std::iter::repeat_n(0, length as usize));
        bytes
    }

    #[test]
    fn test_record_ends_stops_before_a_partial_record() {
        let mut bytes = record(32);
        bytes.extend(record(40));
        bytes.extend(&record(36)[..10]);

        assert_eq!(record_ends(&bytes).unwrap(), vec![36, 80]);
    }

    #[test]
    fn test_record_ends_rejects_a_record_shorter_than_its_fields() {
        assert!(record_ends(&record(8)).is_err());
    }

    #[test]
    fn test_header_length_waits_for_the_whole_header() {
        let mut bytes = BAM_MAGIC.to_vec();
        bytes.extend(3u32.to_le_bytes());
        bytes.extend(b"@HD");
        bytes.extend(1u32.to_le_bytes());
        bytes.extend(5u32.to_le_bytes());
        bytes.extend(b"chr1\0");
        bytes.extend(100u32.to_le_bytes());

        assert_eq!(header_length(&bytes[..20]).unwrap(), None);
        assert_eq!(header_length(&bytes).unwrap(), Some(bytes.len()));
    }

    #[test]
    fn test_header_length_rejects_a_stream_that_is_not_bam() {
        assert!(header_length(b"SAM\x01rest").is_err());
    }
}
