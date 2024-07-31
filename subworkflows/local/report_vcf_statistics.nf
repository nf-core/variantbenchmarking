//
// REPORT_VCF_STATISTICS: SUBWORKFLOW TO REPORT VCF STATS
//

include { SURVIVOR_STATS    } from '../../modules/nf-core/survivor/stats'
include { BCFTOOLS_STATS    } from '../../modules/nf-core/bcftools/stats'

workflow REPORT_VCF_STATISTICS {
    take:
    input_ch    // channel: [val(meta), vcf, index]

    main:

    versions = Channel.empty()

    input_ch.branch{
            def meta = it[0]
            sv:     meta.vartype == "sv" || meta.vartype == "cnv"
            other:  true
        }
        .set{input}

    // use survivor stats to get SV statistics by TYPE
    SURVIVOR_STATS(
        input.sv.map{ meta, vcf, tbi ->
            [ meta, vcf ]
        },
        -1,
        -1,
        -1
    )
    survivor_stats = SURVIVOR_STATS.out.stats
    versions = versions.mix(SURVIVOR_STATS.out.versions.first())

    // use bcftools stats for all files
    BCFTOOLS_STATS(
        input_ch,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    bcftools_stats = BCFTOOLS_STATS.out.stats
    versions = versions.mix(BCFTOOLS_STATS.out.versions.first())

    // Add here a tool, to visualize SV statistics in a histogram.

    emit:
    bcftools_stats
    survivor_stats
    versions
}
