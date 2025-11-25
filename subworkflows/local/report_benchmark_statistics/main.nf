//
// REPORT_BENCHMARK_STATISTICS: SUMMARIZE BENCHMARK REPORTS
//

include { MERGE_REPORTS         } from '../../../modules/local/custom/merge_reports'
include { PLOTS                 } from '../../../modules/local/custom/plots'
include { CREATE_DATAVZRD_INPUT } from '../../../modules/local/custom/create_datavzrd_input'
include { DATAVZRD              } from '../../../modules/nf-core/datavzrd'
include { PLOT_SVLEN_DIST       } from '../../../modules/local/custom/plot_svlen_dist'

workflow REPORT_BENCHMARK_STATISTICS {
    take:
    reports         // channel: [meta, report1, report2, ...]
    evaluations     // channel: [val(meta), vcf.gz, index]
    evaluations_csv // channel: [val(meta), csv]

    main:

    versions = Channel.empty()
    ch_plots = Channel.empty()
    merged_reports = Channel.empty()

    // merge summary statistics from the same benchmarking tool
    MERGE_REPORTS(
        reports
    )
    versions = versions.mix(MERGE_REPORTS.out.versions.first())
    merged_reports = merged_reports.mix(MERGE_REPORTS.out.summary)

    if (!params.skip_plots.contains("metrics")){
        // plot summary statistics
        PLOTS(
            MERGE_REPORTS.out.summary
        )
        ch_plots = ch_plots.mix(PLOTS.out.plots.flatten())
        versions = versions.mix(PLOTS.out.versions.first())
    }

    if (params.variant_type != "snv" && !params.skip_plots.contains("svlenght")){
        // plot INDEL/SV distribution plots
        PLOT_SVLEN_DIST(
            evaluations.groupTuple().mix(evaluations_csv.groupTuple())
        )
        versions = versions.mix(PLOT_SVLEN_DIST.out.versions)
        ch_plots = ch_plots.mix(PLOT_SVLEN_DIST.out.plot)
    }

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
    merged_reports  // channel: [ summary.csv]
}
