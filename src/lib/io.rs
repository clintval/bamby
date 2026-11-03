//! BAM readers and writers over BGZF, single- or multithreaded.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::num::NonZero;
use std::path::Path;

use anyhow::{Context, Result};
use noodles::bam;
use noodles_bgzf as bgzf;
use noodles_bgzf::io::writer::CompressionLevel;

/// A decompressed BAM stream.
pub type Reader = bam::io::Reader<Box<dyn Read>>;

/// Open `path` for reading, decompressing with `threads` workers.
pub fn open(path: &Path, threads: NonZero<usize>) -> Result<Reader> {
    let file = File::open(path).with_context(|| format!("cannot open {}", path.display()))?;
    let inner: Box<dyn Read> = if threads.get() > 1 {
        Box::new(bgzf::io::MultithreadedReader::with_worker_count(
            threads, file,
        ))
    } else {
        Box::new(bgzf::io::Reader::new(BufReader::new(file)))
    };
    Ok(bam::io::Reader::from(inner))
}

/// A BGZF sink that writes the end-of-file block once finished.
pub enum Sink {
    Single(bgzf::io::Writer<Box<dyn Write + Send>>),
    Multi(bgzf::io::MultithreadedWriter<Box<dyn Write + Send>>),
}

impl Sink {
    /// Compress into `path`, or standard output when it is `-`.
    pub fn create(path: &Path, threads: NonZero<usize>, level: CompressionLevel) -> Result<Self> {
        let inner: Box<dyn Write + Send> = if path == Path::new("-") {
            Box::new(BufWriter::new(io::stdout()))
        } else {
            let file =
                File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
            Box::new(BufWriter::new(file))
        };
        Ok(if threads.get() > 1 {
            Self::Multi(
                bgzf::io::multithreaded_writer::Builder::default()
                    .set_compression_level(level)
                    .set_worker_count(threads)
                    .build_from_writer(inner),
            )
        } else {
            Self::Single(
                bgzf::io::writer::Builder::default()
                    .set_compression_level(level)
                    .build_from_writer(inner),
            )
        })
    }

    /// Write the final blocks and flush the underlying writer.
    pub fn finish(self) -> io::Result<()> {
        let mut inner = match self {
            Self::Single(writer) => writer.finish()?,
            Self::Multi(mut writer) => writer.finish()?,
        };
        inner.flush()
    }
}

impl Write for Sink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Single(writer) => writer.write(buf),
            Self::Multi(writer) => writer.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Single(writer) => writer.flush(),
            Self::Multi(writer) => writer.flush(),
        }
    }
}

/// A BAM writer over a [`Sink`].
pub type Writer = bam::io::Writer<Sink>;
