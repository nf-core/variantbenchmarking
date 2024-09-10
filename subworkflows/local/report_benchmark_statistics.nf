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
    report_multiqc  = Channel.empty()

    // merge summary statistics from the same benchmarking tool
    MERGE_REPORTS(
        reports
    )
    versions = versions.mix(MERGE_REPORTS.out.versions.first())
    report_multiqc = report_multiqc.mix(MERGE_REPORTS.out.summary.collect{ meta, summary -> [ summary ] })

    // plot summary statistics
    PLOTS(
        MERGE_REPORTS.out.summary
    )
    versions = versions.mix(PLOTS.out.versions.first())
    //report_multiqc = report_multiqc.mix(PLOTS.out.plots.collect{ meta, plots -> [ plots ] })


    emit:
    versions
    report_multiqc
}
