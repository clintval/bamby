![Banner for bamby](https://raw.githubusercontent.com/clintval/bamby/main/.github/img/bamby_banner.svg)

---

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

## Introduction

The tool `bamby` gives every distinct value of a SAM tag its own read group, so any tool that works per read group or per sample sees one sample per value.
Group single-cell reads by cell barcode, pooled reads by sample barcode, or consensus reads by molecule, all in one BAM instead of one file per value.

## Quick Start

Give every cell barcode in the `CB` tag its own read group:

```bash
bamby -i cells.bam -o grouped.bam -t CB
```

## Per-Cell Variant Calling

`bcftools mpileup` reports each read group's sample as its own column, so every cell is called at once with no per-cell BAMs to split, sort, and merge:

```bash
bamby -i cells.bam -o - -t CB --drop-untagged -l 0 \
    | samtools sort --write-index -o grouped.bam##idx##grouped.bam.bai -

bcftools mpileup -f ref.fa -a AD,DP --per-sample-mF grouped.bam \
    | bcftools call -m -G - -Oz -o cells.vcf.gz
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.
