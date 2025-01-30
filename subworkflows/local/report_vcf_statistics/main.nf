//
// REPORT_VCF_STATISTICS: SUBWORKFLOW TO REPORT VCF STATS
//

include { SURVIVOR_STATS    } from '../../../modules/nf-core/survivor/stats'
include { BCFTOOLS_STATS    } from '../../../modules/nf-core/bcftools/stats'

workflow REPORT_VCF_STATISTICS {
    take:
    input_ch    // channel: [val(meta), vcf, index]

    main:

    versions     = Channel.empty()
    ch_stats     = Channel.empty()

    if (params.variant_type == "structural"){
        // use survivor stats to get SV statistics by TYPE
        SURVIVOR_STATS(
            input_ch.map{ meta, vcf, _tbi ->
                [ meta, vcf ]
            },
            -1,
            -1,
            -1
        )
        ch_stats = ch_stats.mix(SURVIVOR_STATS.out.stats.map{_meta, stats -> stats})
        versions = versions.mix(SURVIVOR_STATS.out.versions.first())
    }


    // use bcftools stats for all files
    BCFTOOLS_STATS(
        input_ch,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]]
    )
    ch_stats = ch_stats.mix(BCFTOOLS_STATS.out.stats.map{_meta, stats -> stats})
    versions = versions.mix(BCFTOOLS_STATS.out.versions.first())

    // TODO: Add here a tool, to visualize SV statistics in a histogram.

    emit:
    ch_stats  // channel: [stats]
    versions  // channel: [versions.yml]
}
