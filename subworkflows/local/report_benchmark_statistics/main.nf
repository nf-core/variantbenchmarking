//
// REPORT_BENCHMARK_STATISTICS: SUMMARIZE BENCHMARK REPORTS
//

include { MERGE_REPORTS         } from '../../../modules/local/custom/merge_reports'
include { PLOTS                 } from '../../../modules/local/custom/plots'
include { CREATE_DATAVZRD_INPUT } from '../../../modules/local/custom/create_datavzrd_input'
include { DATAVZRD              } from '../../../modules/nf-core/datavzrd'

workflow REPORT_BENCHMARK_STATISTICS {
    take:
    reports    // channel: [meta, report1, report2, ...]

    main:

    versions = Channel.empty()
    ch_plots = Channel.empty()

    // merge summary statistics from the same benchmarking tool
    MERGE_REPORTS(
        reports
    )
    versions = versions.mix(MERGE_REPORTS.out.versions.first())

    // plot summary statistics
    PLOTS(
        MERGE_REPORTS.out.summary
    )
    ch_plots = ch_plots.mix(PLOTS.out.plots.flatten())
    versions = versions.mix(PLOTS.out.versions.first())

    MERGE_REPORTS.out.summary
        .map { meta, file -> tuple([vartype: params.variant_type] + [id: meta.benchmark_tool], file) }
        .set { summary }

    // add path to csv file to the datavzrd input
    summary
        .map { meta, summary_file ->
                [ meta, summary_file, file("${projectDir}/assets/datavzrd/${meta.id}.datavzrd.template.yaml", checkIfExists:true) ]
            }
        .set {template}

    CREATE_DATAVZRD_INPUT (
        template
    )

    // use datavzrd to render the report based on the create input
    // input consists of config file and the table itself
    DATAVZRD (
        CREATE_DATAVZRD_INPUT.out.config
    )
    versions = versions.mix(DATAVZRD.out.versions.first())

    emit:
    versions        // channel: [versions.yml]
    ch_plots        // channel: [plots.png]
}
