
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge TP/FP/FN results from different tools.
//

params.options = [:]

include { SURVIVOR_MERGE    } from '../../modules/nf-core/survivor/merge'
include { BCFTOOLS_MERGE    } from '../../modules/nf-core/bcftools/merge'
include { VCF_TO_CSV        } from '../../modules/local/vcf_to_csv'
include { TABIX_BGZIP       } from '../../modules/nf-core/tabix/bgzip'

workflow COMPARE_BENCHMARK_RESULTS {
    take:
    small_ch    // channel: [val(meta), vcf.gz, index]
    sv_ch       // channel: [val(meta), vcf.gz]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:
    versions   = Channel.empty()
    merged_vcfs= Channel.empty()

    // Small Variants
    //
    // MODULE: BCFTOOLS_MERGE
    //
    BCFTOOLS_MERGE(
        small_ch.groupTuple(),
        fasta,
        fai,
        []
    )
    versions = versions.mix(BCFTOOLS_MERGE.out.versions)
    merged_vcfs = merged_vcfs.mix(BCFTOOLS_MERGE.out.merged_variants)

    // SV part
    //
    // MODULE: TABIX_BGZIP
    //
    // unzip vcfs

    TABIX_BGZIP(
        sv_ch
    )
    versions = versions.mix(TABIX_BGZIP.out.versions)

    TABIX_BGZIP.out.output
                .groupTuple()
                .set{vcf_ch}

    //
    // MODULE: SURVIVOR_MERGE
    //
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
    versions = versions.mix(SURVIVOR_MERGE.out.versions)
    merged_vcfs = merged_vcfs.mix(SURVIVOR_MERGE.out.vcf)

    VCF_TO_CSV(
        merged_vcfs,
        fasta,
        fai
    )
    versions = versions.mix(VCF_TO_CSV.out.versions)

    emit:
    versions
}
