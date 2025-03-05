//
// INTERSECT_STATISTICS
//

include { BEDTOOLS_INTERSECT         } from '../../../modules/local/custom/bedtools_intersect'

workflow INTERSECT_STATISTICS {
    take:
    sample_regions
    truth_regions

    main:

    versions        = Channel.empty()

    sample_regions
            .combine(truth_regions)
            .map{test_meta, testvcf, testbed, truthbed -> [test_meta, truthbed, testbed]}
            .set{intersect_ch}

    BEDTOOLS_INTERSECT(
        intersect_ch
    )
    versions      = versions.mix(BEDTOOLS_INTERSECT.out.versions)

    // collect summary reports
    BEDTOOLS_INTERSECT.out.summary
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "intersect"], file) }
        .groupTuple()
        .set{ summary_reports }

    emit:
    versions         // channel: [versions.yml]
    summary_reports  // channel: [meta, summary_report.csv]
}
