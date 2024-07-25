//
// REPORT_VCF_STATISTICS: SUBWORKFLOW TO REPORT VCF STATS
//

include { SURVIVOR_STATS    } from '../../modules/nf-core/survivor/stats'
include { BCFTOOLS_STATS    } from '../../modules/nf-core/bcftools/stats'

workflow REPORT_VCF_STATISTICS {
    take:
    input_ch    // channel: [val(meta), vcf, index]

    main:

    versions=Channel.empty()

    input_ch.branch{
        sv:  it[0].vartype == "sv" || it[0].vartype == "cnv"
        other: true}
        .set{input}

    //
    // SURVIVOR_STATS
    //
    SURVIVOR_STATS(
        input.sv.map{it -> tuple( it[0], it[1])},
        -1,
        -1,
        -1
    )
    survivor_stats = SURVIVOR_STATS.out.stats
    versions = versions.mix(SURVIVOR_STATS.out.versions)

    //
    // BCFTOOLS_STATS
    //
    BCFTOOLS_STATS(
        input_ch,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    bcftools_stats = BCFTOOLS_STATS.out.stats
    versions = versions.mix(BCFTOOLS_STATS.out.versions)

    // Add here a tool, to visualize SV statistics in a histogram.

    emit:
    bcftools_stats
    survivor_stats
    versions
}
