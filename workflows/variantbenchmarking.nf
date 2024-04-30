/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variantbenchmarking_pipeline'

//
// SUBWORKFLOWS: Local Subworkflows
//
include { SV_GERMLINE_BENCHMARK       } from '../subworkflows/local/sv_germline_benchmark'
include { SMALL_GERMLINE_BENCHMARK    } from '../subworkflows/local/small_germline_benchmark'
include { PREPARE_VCFS_TRUTH          } from '../subworkflows/local/prepare_vcfs_truth'
include { PREPARE_VCFS_TEST           } from '../subworkflows/local/prepare_vcfs_test'
include { SV_VCF_CONVERSIONS          } from '../subworkflows/local/sv_vcf_conversion'
include { REPORT_VCF_STATISTICS       } from '../subworkflows/local/report_vcf_statistics'
include { REPORT_BENCHMARK_STATISTICS } from '../subworkflows/local/report_benchmark_statistics'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTBENCHMARKING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_reports       = Channel.empty()
    truth_ch         = Channel.empty()
    high_conf_ch     = Channel.empty()
    bench_ch         = Channel.empty()

    // check mandatory parameters
    println(params.fasta)
    println(params.fai)

    fasta       = Channel.fromPath(params.fasta, checkIfExists: true).map{ it -> tuple([id: it[0].getSimpleName()], it) }.collect()
    fai         = Channel.fromPath(params.fai, checkIfExists: true).map{ it -> tuple([id: it[0].getSimpleName()], it) }.collect()

    // check high confidence files

    truth_small     = params.truth_small        ? Channel.fromPath(params.truth_small, checkIfExists: true).map{ it -> tuple([id: params.sample, vartype:"small"], it) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_small)

    high_conf_small = params.high_conf_small    ? Channel.fromPath(params.high_conf_small, checkIfExists: true).map{ it -> tuple([id: params.sample, vartype:"small"], it) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_small)

    truth_sv        = params.truth_sv           ? Channel.fromPath(params.truth_sv, checkIfExists: true).map{ it -> tuple([id: params.sample, vartype:"sv"], it) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_sv)

    high_conf_sv    = params.high_conf_sv       ? Channel.fromPath(params.high_conf_sv, checkIfExists: true).map{ it -> tuple([id: params.sample, vartype:"sv"], it) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_sv)

    truth_cnv       = params.truth_cnv          ? Channel.fromPath(params.truth_cnv, checkIfExists: true).map{ it -> tuple([id: params.sample, vartype:"cnv"], it) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_cnv)

    high_conf_cnv   = params.high_conf_cnv      ? Channel.fromPath(params.high_conf_cnv, checkIfExists: true).map{ it -> tuple([id: params.sample, vartype:"cnv"], it) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_cnv)

    svync_yaml      = params.svync_yaml         ? Channel.fromPath(params.svync_yaml, checkIfExists: true).collect()
                                                : Channel.empty()


    // TODO: GET FILES FROM IGENOMES ACCORDING TO META.ID

    ch_samplesheet.branch{
            sv:  it[0].vartype == "sv"
            small:  it[0].vartype == "small"
            cnv:  it[0].vartype == "cnv"
            other: false}
            .set{input}

    out_vcf_ch = Channel.empty()

    // PREPROCESSES
    //
    // SUBWORKFLOW: SV_VCF_CONVERSIONS
    //
    // Standardize SV VCFs, tool spesific modifications
    SV_VCF_CONVERSIONS(
        input.sv,
        fasta,
        fai,
        svync_yaml
    )
    ch_versions = ch_versions.mix(SV_VCF_CONVERSIONS.out.versions)
    out_vcf_ch = out_vcf_ch.mix(SV_VCF_CONVERSIONS.out.vcf_ch.map{it -> tuple(it[0], it[1])})
    out_vcf_ch = out_vcf_ch.mix(input.small)
    out_vcf_ch = out_vcf_ch.mix(input.cnv)

    //
    // SUBWORKFLOW: Prepare and normalize input vcfs
    //
    PREPARE_VCFS_TEST(
        out_vcf_ch,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(PREPARE_VCFS_TEST.out.versions)

    PREPARE_VCFS_TRUTH(
        truth_ch,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(PREPARE_VCFS_TRUTH.out.versions)

    // VCF REPORTS AND STATS

    //
    // SUBWORKFLOW: GET STATISTICS OF FILES
    //
    REPORT_VCF_STATISTICS(
        PREPARE_VCFS_TEST.out.vcf_ch.mix(PREPARE_VCFS_TRUTH.out.vcf_ch)
    )
    ch_versions = ch_versions.mix(REPORT_VCF_STATISTICS.out.versions)

    PREPARE_VCFS_TEST.out.vcf_ch.branch{
            sv:  it[0].vartype == "sv"
            small:  it[0].vartype == "small"
            cnv:  it[0].vartype == "cnv"
            other: false}
            .set{test}

    PREPARE_VCFS_TRUTH.out.vcf_ch.branch{
            sv:  it[0].vartype == "sv"
            small:  it[0].vartype == "small"
            cnv:  it[0].vartype == "cnv"
            other: false}
            .set{truth}

    high_conf_ch.branch{
            sv:  it[0].vartype == "sv"
            small:  it[0].vartype == "small"
            cnv:  it[0].vartype == "cnv"
            other: false}
            .set{high_conf}
    // prepare  benchmark set

    if(params.truth_small){
        if(params.high_conf_small){
            test.small.combine(truth.small)
                        .combine(high_conf.small)
                        .map{it -> tuple(it[0], it[1], it[2], it[4], it[5], it[7])}
                        .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.small.combine(truth.small)
                        .map{it -> tuple(it[0], it[1], it[2], it[4], it[5], [])}
                        .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }
    if(params.truth_sv){
        if(params.high_conf_sv){
            test.sv.combine(truth.sv)
                    .combine(high_conf.sv)
                    .map{it -> tuple(it[0], it[1], it[2], it[4], it[5], it[7])}
                    .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.sv.combine(truth.sv)
                    .map{it -> tuple(it[0], it[1], it[2], it[4], it[5], [])}
                    .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }
    if(params.truth_cnv){
        if(params.high_conf_cnv){
            test.cnv.combine(truth.cnv)
                    .combine(high_conf.cnv)
                    .map{it -> tuple(it[0], it[1], it[2], it[4], it[5], it[7])}
                    .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.cnv.combine(truth.cnv)
                    .map{it -> tuple(it[0], it[1], it[2], it[4], it[5], [])}
                    .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }

    // BENCHMARKS
    bench_ch.branch{
            sv:  it[0].vartype == "sv"
            small:  it[0].vartype == "small"
            cnv:  it[0].vartype == "cnv"
            other: false}
            .set{bench_input}
    //
    // SUBWORKFLOW: SMALL_GERMLINE_BENCHMARK
    //
    //Benchmarking spesific to germline samples

    SMALL_GERMLINE_BENCHMARK(
        bench_input.small,
        fasta,
        fai    )
    ch_versions = ch_versions.mix(SMALL_GERMLINE_BENCHMARK.out.versions)
    ch_reports  = ch_reports.mix(SMALL_GERMLINE_BENCHMARK.out.summary_reports)

    //
    // SUBWORKFLOW: SV_GERMLINE_BENCHMARK
    //
    //Benchmarking spesific to germline samples

    SV_GERMLINE_BENCHMARK(
        bench_input.sv,
        fasta,
        fai    )
    ch_versions = ch_versions.mix(SV_GERMLINE_BENCHMARK.out.versions)
    ch_reports  = ch_reports.mix(SV_GERMLINE_BENCHMARK.out.summary_reports)

    // TODO: SOMATIC BENCHMARKING
    if (params.analysis.contains("somatic")){

        // SOMATIC VARIANT BENCHMARKING
        SOMATIC_BENCHMARK(
            bench_input,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(SOMATIC_BENCHMARK.out.versions)
    }

    //
    // SUBWORKFLOW: REPORT_BENCHMARK_STATISTICS
    //
    // Summarize and plot benchmark statistics

    REPORT_BENCHMARK_STATISTICS(
        ch_reports
    )
    ch_versions = ch_versions.mix(REPORT_BENCHMARK_STATISTICS.out.versions)



    // TODO: BENCHMARKING OF CNV
    // https://bioconductor.org/packages/release/bioc/manuals/CNVfilteR/man/CNVfilteR.pdf


    // TODO: TRIO ANALYSIS : MENDELIAN INCONSISTANCE

    // TODO: Compare benchmarking methods!

    //

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
