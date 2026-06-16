//
// REPORT_BENCHMARK_STATISTICS: SUMMARIZE BENCHMARK REPORTS
//

include { MERGE_REPORTS         } from '../../../modules/local/custom/merge_reports'
include { PLOTS_METRICS         } from '../../../modules/local/plots/metrics'
include { DATAVZRD              } from '../../../modules/nf-core/datavzrd'
include { PLOTS_SVLEN_DIST      } from '../../../modules/local/plots/svlen_dist'
include { GAWK as CREATE_DATAVZRD_INPUT } from '../../../modules/nf-core/gawk'

include { REPPY                 } from '../../../modules/local/custom/reppy'

workflow REPORT_BENCHMARK_STATISTICS {
    take:
    reports         // channel: [meta, report1, report2, ...]
    evaluations     // channel: [val(meta), vcf.gz, index]
    evaluations_csv // channel: [val(meta), csv]

    main:

    ch_plots = channel.empty()
    merged_reports = channel.empty()
    evaluations_csv = evaluations_csv.filter { meta, file -> meta.id != "happy" }
    evaluations_csv_happy = evaluations_csv.filter { meta, file -> meta.id == "happy" }

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

    if (params.variant_type == "small" && params.method.contains("happy")){
        // generate HTML report with benchmark statistics using reppy
        // prepare input for reppy
        // you with need a for loop for getting the meta information with the corresponding csv file for each sample and then pass the list of meta and csv file to reppy
        // like a map operator that will take and create a tuple of strings and then concat them into a string.
        evaluations_csv_happy
            .map { meta, file -> "${meta.id}-${meta.caller}:${file}" }
            .collect()
            .map{it -> it.join(" ")}
            .set { reppy_input }
        evaluations_csv_happy
            .map { _meta, file ->  file }
            .collect()
            .set { reppy_input_files}
        REPPY(
            reppy_input_files,
            reppy_input
        )
    }

    MERGE_REPORTS.out.summary
        .map { meta, file -> tuple([vartype: params.variant_type] + [id: meta.benchmark_tool], file) }
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
