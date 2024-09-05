//
// REPORT_BENCHMARK_STATISTICS: SUMMARIZE BENCHMARK REPORTS
//

params.options = [:]

include { MERGE_REPORTS  } from '../../modules/local/merge_reports'          addParams( options: params.options )
include { PLOTS          } from '../../modules/local/plots'                  addParams( options: params.options )
include { CREATE_DATAVZRD_INPUT } from '../../modules/local/create_datavzrd_input'  addParams( options: params.options )
include { DATAVZRD       } from '../../modules/nf-core/datavzrd'           addParams( options: params.options )

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

    // add path to csv file to the datavzrd input
    template = Channel.fromPath( "$projectDir/assets/datavzrd/datavzrd.template.yaml", checkIfExists:true)
    CREATE_DATAVZRD_INPUT (
        template,
        MERGE_REPORTS.out.summary
    )

    // use datavzrd to render the report based on the create input
    // input consists of config file and the table itself
    DATAVZRD (
        CREATE_DATAVZRD_INPUT.out.config
    )
    versions = versions.mix(DATAVZRD.out.versions.first())

    datavzrd_report = DATAVZRD.out.report

    emit:
    versions
    datavzrd_report
}
