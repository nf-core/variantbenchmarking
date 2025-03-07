# nf-core/variantbenchmarking: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.1.0 - [07.03.2025]

Initial release of nf-core/variantbenchmarking, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- CNV benchmarking subworkflow: Truvari (without sequence resolution pctseq = 0) is added as an option.
- _--method intersect_ is implemented enabling intersection two regions (BED) files given. This is especially useful for CNV comparisons where user might only need the segmental matches. The input regions file does not need to be BED file, can also be tool spesfic outputs. According to the tool, formatting will be converted to BED files to be used with bedtools intersect.
- zenodoid added.
- rtgtools vcfeval added for small somatic variant benchmarking with _--squash-ploidy_ parameter.

### `Fixed`

- truth.md links are removed

### `Dependencies`

### `Deprecated`

## 1.0.0 - [24.02.2025]

Initial release of full functioning nf-core/variantbenchmarking

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

## v1.0dev - [12.02.2024]

Initial release of nf-core/variantbenchmarking, created with the [nf-core](https://nf-co.re/) template.
