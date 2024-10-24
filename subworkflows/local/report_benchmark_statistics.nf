//
// REPORT_BENCHMARK_STATISTICS: SUMMARIZE BENCHMARK REPORTS
//

include { MERGE_REPORTS  } from '../../modules/local/merge_reports'
include { PLOTS          } from '../../modules/local/plots'

workflow REPORT_BENCHMARK_STATISTICS {
    take:
    reports    // channel: [meta, report1, report2, ...]

    main:

    versions = Channel.empty()
    // merge summary statistics from the same benchmarking tool
    MERGE_REPORTS(
        reports
    )
    versions = versions.mix(MERGE_REPORTS.out.versions.first())

    // plot summary statistics
    PLOTS(
        MERGE_REPORTS.out.summary
    )
    versions = versions.mix(PLOTS.out.versions.first())


    emit:
    versions
}
