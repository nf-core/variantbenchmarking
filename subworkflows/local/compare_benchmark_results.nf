
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge TP/FP/FN results from different tools.
//

include { SURVIVOR_MERGE    } from '../../modules/nf-core/survivor/merge'
include { BCFTOOLS_MERGE    } from '../../modules/nf-core/bcftools/merge'
include { VCF_TO_CSV        } from '../../modules/local/vcf_to_csv'
include { TABIX_BGZIP       } from '../../modules/nf-core/tabix/bgzip'
include { REFORMAT_HEADER   } from '../../modules/local/reformat_header'

workflow COMPARE_BENCHMARK_RESULTS {
    take:
    small_ch    // channel: [val(meta), vcf.gz, index]
    sv_ch       // channel: [val(meta), vcf.gz, index]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:
    versions    = Channel.empty()
    merged_vcfs = Channel.empty()

    // Small Variants
    REFORMAT_HEADER(
        small_ch
    )
    versions = versions.mix(REFORMAT_HEADER.out.versions.first())

    // merge small variants
    BCFTOOLS_MERGE(
        small_ch.groupTuple(),
        fasta,
        fai,
        []
    )
    versions = versions.mix(BCFTOOLS_MERGE.out.versions.first())
    merged_vcfs = merged_vcfs.mix(BCFTOOLS_MERGE.out.merged_variants)

    // SV part

    // unzip vcfs
    TABIX_BGZIP(
        sv_ch
    )
    versions = versions.mix(TABIX_BGZIP.out.versions.first())

    TABIX_BGZIP.out.output
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

    // convert vcf files to csv
    VCF_TO_CSV(
        merged_vcfs
    )
    versions = versions.mix(VCF_TO_CSV.out.versions.first())

    emit:
    versions
    merged_vcfs
}
