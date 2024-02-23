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
include { SOMATIC_BENCHMARK        } from '../subworkflows/local/somatic_benchmark'
include { SV_GERMLINE_BENCHMARK    } from '../subworkflows/local/sv_germline_benchmark'
include { PREPARE_VCFS_TRUTH       } from '../subworkflows/local/prepare_vcfs_truth'
include { PREPARE_VCFS_TEST        } from '../subworkflows/local/prepare_vcfs_test'
include { SV_VCF_CONVERSIONS       } from '../subworkflows/local/sv_vcf_conversion'
include { REPORT_VCF_STATISTICS as REPORT_STATISTICS_TEST } from '../subworkflows/local/report_vcf_statistics'
include { REPORT_VCF_STATISTICS as REPORT_STATISTICS_TRUTH } from '../subworkflows/local/report_vcf_statistics'

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

    // check mandatory parameters
    println(params.fasta)
    println(params.fai)
    ref         = Channel.fromPath([params.fasta,params.fai], checkIfExists: true).collect()

    // check high confidence files

    truth       = params.truth              ? Channel.fromPath(params.truth, checkIfExists: true).collect()
                                            : Channel.empty()

    high_conf   = params.high_conf          ? Channel.fromPath(params.high_conf, checkIfExists: true).collect()
                                            : Channel.empty()

    svync_yaml  = params.standardization    ? Channel.fromPath("assets/svync/*.yaml").collect()
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
        ref,
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
        ref
    )
    ch_versions = ch_versions.mix(PREPARE_VCFS_TEST.out.versions)

    PREPARE_VCFS_TRUTH(
        truth,
        ref
    )
    ch_versions = ch_versions.mix(PREPARE_VCFS_TRUTH.out.versions)

    // VCF REPORTS AND STATS

    //
    // SUBWORKFLOW: GET STATISTICS OF FILES
    //
    //REPORT_STATISTICS_TRUTH(
    //
    REPORT_STATISTICS_TEST(
        PREPARE_VCFS_TEST.out.vcf_ch
    )
    REPORT_STATISTICS_TRUTH(
        PREPARE_VCFS_TRUTH.out.vcf_ch
    )
    ch_versions = ch_versions.mix(PREPARE_VCFS_TRUTH.out.versions)

    // prepare  benchmark set
    if (params.high_conf){
        PREPARE_VCFS_TEST.out.vcf_ch.combine(PREPARE_VCFS_TRUTH.out.vcf_ch.map { it -> tuple(it[1], it[2]) })
                                .combine(high_conf)
                                .set{bench_ch}
    }else{
        PREPARE_VCFS_TEST.out.vcf_ch.combine(PREPARE_VCFS_TRUTH.out.vcf_ch.map { it -> tuple(it[1], it[2]) })
                                .map{it -> tuple(it[0], it[1], it[2],it[3],it[4],[])}
                                .set{bench_ch}
    }

    // BENCHMARKS

    bench_ch.branch{
            sv:  it[0].vartype == "sv"
            small:  it[0].vartype == "small"
            cnv:  it[0].vartype == "cnv"
            other: false}
            .set{input}

    //
    // SUBWORKFLOW: SV_GERMLINE_BENCHMARK
    //
    //Benchmarking spesific to germline samples

    SV_GERMLINE_BENCHMARK(
        input.sv,
        ref    )
    ch_versions = ch_versions.mix(SV_GERMLINE_BENCHMARK.out.versions)

    // TODO: SOMATIC EBNCHMARKING
    if (params.analysis.contains("somatic")){

        // SOMATIC VARIANT BENCHMARKING
        SOMATIC_BENCHMARK(
            bench_ch,
            ref
        )
        ch_versions = ch_versions.mix(SOMATIC_BENCHMARK.out.versions)
    }

    // TODO: NEED A TOOL TO COLLECT METRICS AND ROCS LIKE DATAVZRD OR SQLITE DATABASE


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
