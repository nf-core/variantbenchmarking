//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

params.options = [:]

include { SURVIVOR_STATS    } from '../../modules/nf-core/survivor/stats'      addParams( options: params.options )

workflow REPORT_STATISTICS {
    take:
    input_ch    // channel: [val(meta), vcf, index]

    main:

    versions=Channel.empty()

    //
    // SURVIVOR_STATS
    //

    SURVIVOR_STATS(
        input_ch,
        -1,
        -1,
        -1
    )
    stats = SURVIVOR_STATS.out.stats
    versions = versions.mix(SURVIVOR_STATS.out.versions)


    emit:
    stats
    versions
}
