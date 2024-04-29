//
// REPORT_BENCHMARK_STATISTICS: SUMMARIZE BENCHMARK REPORTS
//

params.options = [:]

include { MERGE_REPORTS    } from '../../modules/local/merge_reports'       addParams( options: params.options )

workflow REPORT_BENCHMARK_STATISTICS {
    take:
    reports    // channel: [val(benchmark_tool), report1, report2, ...]

    main:

    versions=Channel.empty()

    MERGE_REPORTS(reports)

    versions = versions.mix(MERGE_REPORTS.out.versions)

    emit:
    versions
}
