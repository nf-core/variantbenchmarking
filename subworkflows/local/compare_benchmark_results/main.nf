
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge TP/FP/FN results from different tools.
//

include { SURVIVOR_MERGE    } from '../../../modules/nf-core/survivor/merge'
include { BCFTOOLS_MERGE    } from '../../../modules/nf-core/bcftools/merge'
include { VCF_TO_CSV        } from '../../../modules/local/custom/vcf_to_csv'
include { REFORMAT_HEADER   } from '../../../modules/local/custom/reformat_header'
include { TABIX_BGZIP as TABIX_BGZIP_UNZIP } from '../../../modules/nf-core/tabix/bgzip'

workflow COMPARE_BENCHMARK_RESULTS {
    take:
    evaluations // channel: [val(meta), vcf.gz, index]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:
    versions    = Channel.empty()
    merged_vcfs = Channel.empty()

    if (params.variant_type == "small" | params.variant_type == "snv" | params.variant_type == "indel"){

        // Small Variants
        REFORMAT_HEADER(
            evaluations
        )
        versions = versions.mix(REFORMAT_HEADER.out.versions.first())

        // merge small variants
        BCFTOOLS_MERGE(
            evaluations.groupTuple(),
            fasta,
            fai,
            []
        )
        versions = versions.mix(BCFTOOLS_MERGE.out.versions.first())
        merged_vcfs = merged_vcfs.mix(BCFTOOLS_MERGE.out.merged_variants)
    }
    else{
        // SV part
        // unzip vcfs
        TABIX_BGZIP_UNZIP(
            evaluations
        )
        versions = versions.mix(TABIX_BGZIP_UNZIP.out.versions.first())

        TABIX_BGZIP_UNZIP.out.output
            .groupTuple()
            .set{vcf_ch}

        // Merge Benchmark SVs from different tools
        SURVIVOR_MERGE(
            vcf_ch,
            1000,
            1,
            1,
            0,
            0,
            30
        )
        versions = versions.mix(SURVIVOR_MERGE.out.versions.first())
        merged_vcfs = merged_vcfs.mix(SURVIVOR_MERGE.out.vcf)
    }

    // convert vcf files to csv
    VCF_TO_CSV(
        merged_vcfs
    )
    versions = versions.mix(VCF_TO_CSV.out.versions.first())

    emit:
    merged_vcfs  // channel: [val(meta), vcf]
    versions     // channel: [versions.yml]
}
