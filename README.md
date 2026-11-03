# bamby

[![Install with bioconda](https://img.shields.io/badge/Install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/bamby/README.html)
[![Anaconda Version](https://anaconda.org/bioconda/bamby/badges/version.svg)](http://bioconda.github.io/recipes/bamby/README.html)
[![Build Status](https://github.com/clintval/bamby/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/bamby/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/clintval/bamby/badge.svg?branch=main)](https://coveralls.io/github/clintval/bamby?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Group BAM records into read groups by the value of a tag.

Install with mamba, conda, or run directly with pixi:

```bash
pixi exec \
    -c conda-forge -c bioconda \
    bamby --help
```

## Quick Start

Give every cell barcode in `CB` its own read group:

```bash
bamby -i cells.bam -o grouped.bam
```

Every tool that works per read group or per sample now sees one sample per cell, with no need to split the BAM into one file per cell.
For example, `bcftools mpileup` reports each cell as its own sample column:

```bash
bamby -i cells.bam -o - -t CB --drop-untagged -l 0 \
| samtools sort --write-index -o grouped.bam##idx##grouped.bam.bai

bcftools mpileup -f ref.fa -a AD,DP grouped.bam | bcftools call -m -G - -Oz -o cells.vcf.gz
```

## How It Works

The input is read twice: once to learn every value of the tag and build the header, then again to write each record with its new `RG` tag.
Records keep their order, so a coordinate-sorted input stays coordinate-sorted.

Each new read group is named for the record's read group and the tag value:

| Field | Value |
| --- | --- |
| `ID` | The record's read group and the value, joined by a period (`A.AACGTT`), or the value alone when the record has no read group. |
| `SM` | The value. |
| Others | Every other field (`LB`, `PL`, `PU`, ...) copied from the record's read group. |

Records without the tag keep their read group, or are dropped with `--drop-untagged`.
A read group that no written record uses any more is removed from the header, and new read groups are listed in sorted order.
String, character, hex, and integer tags can name read groups; float and array tags cannot.

## Options

| Option | Default | Description |
| --- | --- | --- |
| `-i`, `--input` | | Input BAM, read twice, so it must be a file. |
| `-o`, `--output` | | Output BAM, or `-` for standard output. |
| `-t`, `--tag` | `CB` | The tag whose values name the read groups. |
| `-d`, `--drop-untagged` | | Drop records without the tag. |
| `-@`, `--threads` | `1` | BGZF worker threads, for reading and again for writing. |
| `-l`, `--compression-level` | `6` | BGZF compression level of the output, from 0 (none) to 12 (best). |

## Related Tools

- [`samtools split -d TAG`](https://www.htslib.org/doc/samtools-split.html) writes one file per tag value.
- [`sinto tagtorg`](https://timoast.github.io/sinto/) does the same job in Python.
