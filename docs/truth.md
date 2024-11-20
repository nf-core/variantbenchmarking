# nf-core/variantbenchmarking: Truth files

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/variantbenchmarking/truth](https://nf-co.re/variantbenchmarking/truth)

## Defining Truth VCF and High confidence BED files

This pipeline requires a set of Truth VCF, as a baseline for comparisons, and a high confidence bed files, to restrict analysis to regions. Although, those sets can be anything depending on the type of the analysis, for benchmarking of human genomes there are golden set of samples provided by [Genome in a Bottle project](https://www.nist.gov/programs-projects/genome-bottle) and [SEQC2 consortium](https://sites.google.com/view/seqc2/home/data-analysis/high-confidence-somatic-snv-and-indel-v1-2).

Below, please find some set truths can be used for analysis:

- params.analysis == 'germline'

  - params.genome == 'GRCh38'

    - params.vartype == 'small'

                    truth_id          = "HG002"
                    truth_vcf         = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
                    regions_bed       = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

  - params.genome == 'GRCh37'

    - params.vartype == 'small'

                    truth_id          = "HG002"
                    truth_vcf         = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/SupplementaryFiles/HG002_GRCh37_1_22_v4.2.1_highconf.vcf.gz"
                    regions_bed       = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/SupplementaryFiles/HG002_GRCh37_1_22_v4.2.1_highconf.bed"

    - params.vartype == 'structural'

                    truth_id          = "HG002"
                    truth_vcf          = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz"
                    regions_bed        = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/HG002_SVs_Tier1_v0.6.bed"

- params.analysis == 'somatic'

  - params.genome == 'GRCh38'

    - params.vartype == 'structural'

                    truth_id          = "HG002"
                    truth_vcf         = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.vcf.gz"
                    regions_bed       = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/GRCh38/StructuralVariant/HG002_GRCh38_CMRG_SV_v1.00.bed"

    - params.vartype == 'snv'

                    truth_id          = "SEQC2"
                    truth_vcf         = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz"
                    regions_bed       = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/High-Confidence_Regions_v1.2.bed"

    - params.vartype == 'indel'

                    truth_id          = "SEQC2"
                    truth_vcf         = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz"
                    regions_bed       = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/High-Confidence_Regions_v1.2.bed"
