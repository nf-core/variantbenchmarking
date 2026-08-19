//
// REPORT_BENCHMARK_STATISTICS: SUMMARIZE BENCHMARK REPORTS
//

include { MERGE_REPORTS         } from '../../../modules/local/custom/merge_reports'
include { PLOTS_METRICS         } from '../../../modules/local/plots/metrics'
include { DATAVZRD              } from '../../../modules/nf-core/datavzrd'
include { PLOTS_SVLEN_DIST      } from '../../../modules/local/plots/svlen_dist'
include { GAWK as CREATE_DATAVZRD_INPUT } from '../../../modules/nf-core/gawk'

workflow REPORT_BENCHMARK_STATISTICS {
    take:
    reports         // channel: [meta, report1, report2, ...]
    evaluations     // channel: [val(meta), vcf.gz, index]
    evaluations_csv // channel: [val(meta), csv]

    main:

    ch_plots = channel.empty()
    merged_reports = channel.empty()

    // merge summary statistics from the same benchmarking tool
    MERGE_REPORTS(
        reports
    )
    merged_reports = merged_reports.mix(MERGE_REPORTS.out.summary)

    if (!params.skip_plots.contains("metrics")){
        // plot summary statistics
        PLOTS_METRICS(
            MERGE_REPORTS.out.summary
        )
        ch_plots = ch_plots.mix(PLOTS_METRICS.out.plots.flatten())
    }

    if (params.variant_type != "snv" && !params.skip_plots.contains("svlength")){
        // plot INDEL/SV distribution plots
        evaluations.map { item ->
            tuple(item[0], item[1])
        }.groupTuple()
        .mix(
            evaluations_csv.map { item ->
                tuple(item[0], item[1])
            }.groupTuple()
        )
        .set { svlen_input }

        PLOTS_SVLEN_DIST(svlen_input)
        ch_plots = ch_plots.mix(PLOTS_SVLEN_DIST.out.plot)
    }

    MERGE_REPORTS.out.summary
        .map { meta, report -> tuple([vartype: params.variant_type] + [id: meta.benchmark_tool], report) }
        .set { summary }

    // add path to csv file to the datavzrd input
    summary
        .map { meta, summary_file ->
            def updated_meta = meta + [ csv: summary_file.toString() ]
            def template_file = file("${projectDir}/assets/datavzrd/${meta.id}.datavzrd.template.yaml", checkIfExists: true)
            [ updated_meta, template_file ]
        }
        .set { template_ch }

    CREATE_DATAVZRD_INPUT (
        template_ch,
        [],
        false
    )

    CREATE_DATAVZRD_INPUT.out.output
        .map { meta, yaml_file ->
            def clean_meta = meta.findAll { it.key != 'csv' }
            [ clean_meta, yaml_file ]
        }
        .set { clean_datavzrd_input_ch }

    DATAVZRD (
        clean_datavzrd_input_ch.join(summary)
    )

    emit:
    ch_plots        // channel: [ plots.png ]
    merged_reports  // channel: [ meta, summary.csv]
}
