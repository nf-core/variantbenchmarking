# nf-core/variantbenchmarking: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/variantbenchmarking/usage](https://nf-co.re/variantbenchmarking/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the test vcf you would like to analyze before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 4 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Full samplesheet

The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

```csv title="samplesheet.csv"
id,test_vcf,caller,vartype
test1,test1.vcf.gz,delly,sv
test2,test2.vcf,gatk,small
test3,test3.vcf.gz,cnvkit,cnv
```

| Column     | Description                                                                                                                                                                            |
| ---------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `id`       | Custom id name per test vcf. This entry will be identical.                                                                                                                             |
| `test_vcf` | The VCF file to use as benchmarking test input. The same file can be used in more than one row. File can be either vcf or vcf.gz.                                                      |
| `caller`   | Variant caller method used to generate test VCF file. There can be more than one test vcf for the same caller. For unknown caller use 'unknown'                                        |
| `vartype`  | Variant type to apply benchmarking. Variant type can be only one of these: small, sv, snv, indel and cnv.                                                                              |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.


## Defining Truth VCF and High confidence BED files

The following parameters has to be defined for each type of benchmarking analysis. The following parameters defined the exact paths to the truth files:

- `--sample`: Sample parameter defines the same name of the truth set. Examples: `HG002`, `SEQC2`, `HG001`, `HG003`, `CHM13`.
- `--analysis`: The type of analysis to perform: `germline` or `somatic`.
- `--method`: List the benchmarking methods to apply. By default all available tools will be applied according to the variant types provided. Available tools: `truvari, svanalyzer, happy, sompy, rtgtools, wittyer`.

*Small variant benchmarking:*

- `--truth_small`: Path to the golden set VCF files combined for SNVs and indels, required for germline benchmarking (vcf or vcf.gz)
- `--high_conf_small`: Path to the high confidence BED files for SNVs and indels, required for germline benchmarking  (bed or bed.gz)
- `--truth_snv`: Path to the golden set VCF files for SNVs, required for somatic benchmarking (vcf or vcf.gz)
- `--high_conf_snv`: Path to the high confidence BED files for SNVs, required for somatic benchmarking (bed or bed.gz)
- `--truth_indel`: Path to the golden set VCF files for indels, required for somatic benchmarking (vcf or vcf.gz)
- `--high_conf_indel`: Path to the high confidence BED files for indels, required for somatic benchmarking (bed or bed.gz)

*Structural variant benchmarking:*

- `--truth_sv`: Path to the golden set VCF files for SVs, required for germline and somatic benchmarking (vcf or vcf.gz)
- `--high_conf_sv`: Path to the high confidence BED files for SVs, required for germline and somatic benchmarking (bed or bed.gz)

*Copy Number Variation benchmarking:*

- `--truth_cnv`: Path to the golden set VCF files for CNVs, required for germline and somatic benchmarking (vcf or vcf.gz)
- `--high_conf_cnv`: Path to the high confidence BED files for CNVs, required for germline and somatic benchmarking (bed or bed.gz)

*Using truth.config*

`conf/truth.config` file contains some readily available truth files for germline and somatic analysis. In order to activate usage one has to

1. use `--genome` [`GRCh37` or `GRCh38`]
2. define `--sample` [`HG002` or `SEQC2`]
3. turn off `--itruth_ignore false`

## Lifting over truth sets

This workflow comes with a liftover option for truth sets. In order to activate liftover use `--liftover true`.

- `--chain`: This workflow uses picard tools for lifting over and a chain file has to be provided specific to the input truth vcf. Some examples can be found [here](https://genome.ucsc.edu/goldenPath/help/chain.html)
- `--rename_chr`: Renaming chromosomes is required after liftover process. Some examples can be found under `assets/rename_contigs` directory.

Note: these two files are also provided under `itruth.config`. An example usage can be found in `conf/test_liftover.config`

## Standardization and normalization parameters

Consistent formatting and alignment of variants in test and truth VCF files for accurate comparison is controlled by *sv_standardization* and *preprocesses*.

- `--sv_standardization`: The standardization methods to perform on the input files. Should be a comma-separated list of one or more of the following options: `homogenize,svync`.
  - `homogenize`: makes use of [variant-extractor](https://github.com/EUCANCan/variant-extractor)
  - `svync`: makes use of [svync](https://github.com/nvnieuwk/svync)

- `--preprocesses`: The preprocessing steps to perform on the input files. Should be a comma-separated list of one or more of the following options: `normalization,deduplication,prepy,filter_contigs`
  - `normalization`: Splits multi-allelic variants in test and truth VCF files ([bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm))
  - `deduplication`: Deduplicates variants in test and truth VCF files ([bcftools norm](https://samtools.github.io/bcftools/bcftools.html#norm))
  - `prepy`: Uses prepy in order to normalize test files. This option is only applicable for happy benchmarking of germline analysis ([prepy](https://github.com/Illumina/hap.py/tree/master))
  - `filter_contigs`: Filter out extra contigs. It is common for truth files not to include extra contigs.

## Using multi-sample vcf inputs

If the input test vcf contains more than one sample, then user has to define which sample name to use. `subsample` will added to the samplesheet as an additional column as follows:


```csv title="samplesheet.csv"
id,test_vcf,caller,vartype,subsample
test1,test1.vcf.gz,delly,sv,"TUMOR"
test2,test2.vcf,gatk,small,"NA128120"
test3,test3.vcf.gz,cnvkit,cnv,
```

Note that, this option can be inevitable for somatic analysis since most of the callers reports both normal and tumor genotypes in the same vcf file.

## Filtering parameters

- `--exclude_expression`: Use [bcftools expressions](https://samtools.github.io/bcftools/bcftools.html#expressions) to exclude variants. Default:null
- `--include_expression`: Use [bcftools expressions](https://samtools.github.io/bcftools/bcftools.html#expressions) to include variants. Default:null

*Parameters applicable only to Structural Variants*

- `--min_sv_size`: Minimum SV size of variants to benchmark, 0 to disable , Default:30
- `--max_sv_size`: Maximum SV size of variants to benchmark, -1 to disable , Default:-1
- `--min_allele_freq`: Minimum Alele Frequency of variants to benchmark, Use -1 to disable , Default:-1
- `--min_num_reads`: Minimum number of read supporting variants to benchmark, Use, -1 to disable , Default:-1

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/variantbenchmarking --input ./samplesheet.csv --outdir ./results -profile docker --genome GRCh37 --sample HG002 --analysis germline
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/variantbenchmarking -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
sample: 'HG002'
analysis: 'germline'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/variantbenchmarking
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/variantbenchmarking releases page](https://github.com/nf-core/variantbenchmarking/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `test_liftover`
  - A profile with a complete configuration for using liftover of HG002 hg38 truth set to hg37
  - Includes links to test data so needs no other parameters
- `test_germline`
  - A profile with a complete configuration for a full test of HG002 sample from germline analysis
  - Includes links to test data so needs no other parameters
- `test_somatic`
  - A profile with a complete configuration for a full test of SEQC2 sample from somatic analysis
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
