# nf-core/variantbenchmarking: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 1.3.0 - dev

### `Added`

- Add Precision vs recall plot to compare benchmarking tools ([#198](https://github.com/nf-core/variantbenchmarking/pull/198))
- Adding a test profile to showcase running with GA4GH best practices happy with rtgtools engine([#189](https://github.com/nf-core/variantbenchmarking/pull/189))
- nf-core-template-merge-3.2.1 ([#193](https://github.com/nf-core/variantbenchmarking/pull/193))
- Add `rtg bndeval` for SVTPE=BND benchmarking ([#195](https://github.com/nf-core/variantbenchmarking/pull/195))
- Split benchmarking VCF for happy properly to create comparison tables ([#196](https://github.com/nf-core/variantbenchmarking/pull/196))
- Documentation about test configs in conf/tests ([#203](https://github.com/nf-core/variantbenchmarking/pull/203))
- Template update to v3.3.1 [#200](https://github.com/nf-core/variantbenchmarking/pull/200)
- Template update to v3.3.2 [#205](https://github.com/nf-core/variantbenchmarking/pull/205)

### `Fixed`

- Fixing bedtools intersection logic for CNV calculations ([#192](https://github.com/nf-core/variantbenchmarking/pull/192))
- Seperate subworkflows for each benchmark method ([#203](https://github.com/nf-core/variantbenchmarking/pull/203))
- Making test cases more realistic and adding docs how to use different tools in different context ([#204](https://github.com/nf-core/variantbenchmarking/pull/204))
- #205 overwrites #192 making nf-schema to latest, update metromap [#206](https://github.com/nf-core/variantbenchmarking/pull/206)

### `Dependencies`

- Downgrade nf-schema to fix CI tests ([#192](https://github.com/nf-core/variantbenchmarking/pull/192))

### `Deprecated`

## 1.2.0 - [31.03.2025]

### `Added`

- Updated metromap
- [Implementing optional bed files for benchmarking tools (happy, rtgtools, sompy) & adding stratifications for genome](https://github.com/nf-core/variantbenchmarking/pull/167)
  - added regions_bed vs targets_bed. regions_bed works as the same way Bcftools -R and targets_bed as Bcftools -T. This differenciation only exist for happy, rtgtools and sompy. truvari, svanalyzer and wittyer uses only regions_bed.

### `Fixed`

- [official logo cleanup](https://github.com/nf-core/variantbenchmarking/pull/171)
- [bcftools_reheader issue for HPCs](https://github.com/nf-core/variantbenchmarking/pull/170)
- [json-schema-improvements](https://github.com/nf-core/variantbenchmarking/pull/168)

### `Dependencies`

### `Deprecated`

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
