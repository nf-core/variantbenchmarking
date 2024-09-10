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
include { SUBSAMPLE_VCF_TEST          } from '../subworkflows/local/subsample_vcf_test'
include { PREPARE_VCFS_TRUTH          } from '../subworkflows/local/prepare_vcfs_truth'
include { PREPARE_VCFS_TEST           } from '../subworkflows/local/prepare_vcfs_test'
include { SV_VCF_CONVERSIONS          } from '../subworkflows/local/sv_vcf_conversion'
include { REPORT_VCF_STATISTICS       } from '../subworkflows/local/report_vcf_statistics'
include { SV_GERMLINE_BENCHMARK       } from '../subworkflows/local/sv_germline_benchmark'
include { SMALL_GERMLINE_BENCHMARK    } from '../subworkflows/local/small_germline_benchmark'
include { CNV_GERMLINE_BENCHMARK      } from '../subworkflows/local/cnv_germline_benchmark'
include { SMALL_SOMATIC_BENCHMARK     } from '../subworkflows/local/small_somatic_benchmark'
include { REPORT_BENCHMARK_STATISTICS } from '../subworkflows/local/report_benchmark_statistics'
include { COMPARE_BENCHMARK_RESULTS   } from '../subworkflows/local/compare_benchmark_results'

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
    ch_summary       = Channel.empty()
    truth_ch         = Channel.empty()
    high_conf_ch     = Channel.empty()
    bench_ch         = Channel.empty()
    sv_evals_ch      = Channel.empty()
    small_evals_ch   = Channel.empty()

    //// create reference channels ////

    fasta       = Channel.fromPath(params.fasta, checkIfExists: true)
                    .map{ fasta -> tuple([id: fasta.getSimpleName()], fasta) }.collect()
    fai         = Channel.fromPath(params.fai, checkIfExists: true)
                    .map{ fai -> tuple([id: fai.getSimpleName()], fai) }.collect()

    //// check high confidence files ////

    // Germline

    truth_small     = params.truth_small        ? Channel.fromPath(params.truth_small, checkIfExists: true).map{ vcf -> tuple([id: params.sample, vartype:"small"], vcf) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_small)

    high_conf_small = params.high_conf_small    ? Channel.fromPath(params.high_conf_small, checkIfExists: true).map{ bed -> tuple([id: params.sample, vartype:"small"], bed) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_small)

    truth_sv        = params.truth_sv           ? Channel.fromPath(params.truth_sv, checkIfExists: true).map{ vcf -> tuple([id: params.sample, vartype:"sv"], vcf) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_sv)

    high_conf_sv    = params.high_conf_sv       ? Channel.fromPath(params.high_conf_sv, checkIfExists: true).map{ bed -> tuple([id: params.sample, vartype:"sv"], bed) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_sv)

    truth_cnv       = params.truth_cnv          ? Channel.fromPath(params.truth_cnv, checkIfExists: true).map{ vcf -> tuple([id: params.sample, vartype:"cnv"], vcf) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_cnv)

    high_conf_cnv   = params.high_conf_cnv      ? Channel.fromPath(params.high_conf_cnv, checkIfExists: true).map{ bed -> tuple([id: params.sample, vartype:"cnv"], bed) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_cnv)

    // Somatic
    // snv and indel seperation only possible for somatic cases

    truth_snv       = params.truth_snv          ? Channel.fromPath(params.truth_snv, checkIfExists: true).map{ vcf -> tuple([id: params.sample, vartype:"snv"], vcf) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_snv)

    high_conf_snv   = params.high_conf_snv      ? Channel.fromPath(params.high_conf_snv, checkIfExists: true).map{ bed -> tuple([id: params.sample, vartype:"snv"], bed) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_snv)

    truth_indel     = params.truth_indel        ? Channel.fromPath(params.truth_indel, checkIfExists: true).map{ vcf -> tuple([id: params.sample, vartype:"indel"], vcf) }.collect()
                                                : Channel.empty()
    truth_ch        = truth_ch.mix(truth_indel)

    high_conf_indel = params.high_conf_indel    ? Channel.fromPath(params.high_conf_indel, checkIfExists: true).map{ bed -> tuple([id: params.sample, vartype:"indel"], bed) }.collect()
                                                : Channel.empty()
    high_conf_ch    = high_conf_ch.mix(high_conf_indel)


    // SDF file for RTG-tools eval
    sdf             = params.sdf                ? Channel.fromPath(params.sdf, checkIfExists: true).map{ sdf -> tuple([id: sdf.getSimpleName()], sdf) }.collect()
                                                : Channel.empty()

    // read chainfile, liftover genome and rename chr files if liftover is true
    chain           = Channel.empty()
    rename_chr      = Channel.empty()
    dictionary      = Channel.empty()

    if (params.liftover){
        chain           = params.chain          ? Channel.fromPath(params.chain, checkIfExists: true).map{ bed -> tuple([id: bed.getSimpleName()], bed) }.collect()
                                                : Channel.empty()

        rename_chr      = params.rename_chr     ? Channel.fromPath(params.rename_chr, checkIfExists: true).map{ txt -> tuple([id: txt.getSimpleName()], txt) }.collect()
                                                : Channel.empty()

        dictionary      = params.dictionary     ? Channel.fromPath(params.dictionary, checkIfExists: true).map{ dict -> tuple([id: dict.getSimpleName()], dict) }.collect()
                                                : Channel.empty()
    }
    // PREPROCESSES

    // subsample multisample vcf if necessary
    ch_samplesheet.branch{
            def meta = it[0]
            multisample: meta.subsample != null
            other: true}
        .set{input}

    out_vcf_ch  = Channel.empty()

    SUBSAMPLE_VCF_TEST(
        input.multisample
    )
    ch_versions = ch_versions.mix(SUBSAMPLE_VCF_TEST.out.versions)
    out_vcf_ch  = out_vcf_ch.mix(SUBSAMPLE_VCF_TEST.out.vcf_ch,
                                input.other)
    vcf_ch      = out_vcf_ch

    // Branch out according to the analysis
    vcf_ch.branch{
            sv:  it[0].vartype == "sv"
            other: true}
            .set{test_ch}

    out_vcf_ch = Channel.empty()

    // Standardize SV VCFs, tool spesific modifications
    SV_VCF_CONVERSIONS(
        test_ch.sv,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(SV_VCF_CONVERSIONS.out.versions)
    out_vcf_ch = out_vcf_ch.mix(SV_VCF_CONVERSIONS.out.vcf_ch.map{it -> tuple(it[0], it[1])},
                            test_ch.other)

    // Prepare and normalize input vcfs
    PREPARE_VCFS_TEST(
        out_vcf_ch,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(PREPARE_VCFS_TEST.out.versions)

    // Prepare and normalize truth vcfs
    PREPARE_VCFS_TRUTH(
        truth_ch,
        high_conf_ch,
        fasta,
        fai,
        chain,
        rename_chr,
        dictionary
    )
    high_conf_ch = PREPARE_VCFS_TRUTH.out.high_conf_ch
    ch_versions = ch_versions.mix(PREPARE_VCFS_TRUTH.out.versions)

    // VCF REPORTS AND STATS

    // get statistics for normalized input files
    REPORT_VCF_STATISTICS(
        PREPARE_VCFS_TEST.out.vcf_ch.mix(PREPARE_VCFS_TRUTH.out.vcf_ch)
    )
    ch_versions = ch_versions.mix(REPORT_VCF_STATISTICS.out.versions)

    // branch out input test files
    PREPARE_VCFS_TEST.out.vcf_ch.branch{
            def meta = it[0]
            sv:     meta.vartype == "sv"
            small:  meta.vartype == "small"
            cnv:    meta.vartype == "cnv"
            snv:    meta.vartype == "snv"
            indel:  meta.vartype == "indel"
            other:  false
        }
        .set{test}

    // branch out truth vcf files
    PREPARE_VCFS_TRUTH.out.vcf_ch.branch{
            def meta = it[0]
            sv:     meta.vartype == "sv"
            small:  meta.vartype == "small"
            cnv:    meta.vartype == "cnv"
            snv:    meta.vartype == "snv"
            indel:  meta.vartype == "indel"
            other:  false
        }
        .set{truth}

    // branch out high confidence bed files
    high_conf_ch.branch{
            def meta = it[0]
            sv:     meta.vartype == "sv"
            small:  meta.vartype == "small"
            cnv:    meta.vartype == "cnv"
            snv:    meta.vartype == "snv"
            indel:  meta.vartype == "indel"
            other:  false
        }
        .set{high_conf}

    // prepare  benchmark sets
    if(params.truth_small){
        if(params.high_conf_small){
            test.small.combine(truth.small)
                .combine(high_conf.small)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi, high_meta, high_bed ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, high_bed ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.small.combine(truth.small)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, [] ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }
    if(params.truth_sv){
        if(params.high_conf_sv){
            test.sv.combine(truth.sv)
                .combine(high_conf.sv)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi, high_meta, high_bed ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, high_bed ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.sv.combine(truth.sv)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, [] ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }
    if(params.truth_cnv){
        if(params.high_conf_cnv){
            test.cnv.combine(truth.cnv)
                .combine(high_conf.cnv)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi, high_meta, high_bed ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, high_bed ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.cnv.combine(truth.cnv)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, [] ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }
    if(params.truth_snv){
        if(params.high_conf_snv){
            test.snv.combine(truth.snv)
                .combine(high_conf.snv)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi, high_meta, high_bed ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, high_bed ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.snv.combine(truth.snv)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, [] ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }
    if(params.truth_indel){
        if(params.high_conf_indel){
            test.indel.combine(truth.indel)
                .combine(high_conf.indel)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi, high_meta, high_bed ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, high_bed ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
        else{
            test.indel.combine(truth.indel)
                .map{ test_meta, test_vcf, test_tbi, truth_meta, truth_vcf, truth_tbi ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, [] ]
                }
                .set{bench}
            bench_ch = bench_ch.mix(bench)
        }
    }

    // branch out combined benchmark sets
    bench_ch.branch{
            def meta = it[0]
            sv:     meta.vartype == "sv"
            small:  meta.vartype == "small"
            cnv:    meta.vartype == "cnv"
            snv:    meta.vartype == "snv"
            indel:  meta.vartype == "indel"
            other:  false
        }
        .set{bench_input}


    // Perform SV benchmarking - for now it also works for somatic cases!
    SV_GERMLINE_BENCHMARK(
        bench_input.sv,
        fasta,
        fai
    )
    ch_versions = ch_versions.mix(SV_GERMLINE_BENCHMARK.out.versions)
    ch_reports  = ch_reports.mix(SV_GERMLINE_BENCHMARK.out.summary_reports)
    sv_evals_ch = sv_evals_ch.mix(SV_GERMLINE_BENCHMARK.out.tagged_variants)
    ch_summary  = ch_summary.mix(SV_GERMLINE_BENCHMARK.out.report_multiqc)

    if (params.analysis.contains("germline")){

        // Benchmarking spesific to small germline samples
        SMALL_GERMLINE_BENCHMARK(
            bench_input.small,
            fasta,
            fai,
            sdf
        )
        ch_versions = ch_versions.mix(SMALL_GERMLINE_BENCHMARK.out.versions)
        ch_reports  = ch_reports.mix(SMALL_GERMLINE_BENCHMARK.out.summary_reports)
        small_evals_ch = small_evals_ch.mix(SMALL_GERMLINE_BENCHMARK.out.tagged_variants)
        ch_summary  = ch_summary.mix(SMALL_GERMLINE_BENCHMARK.out.report_multiqc)

        // Benchmarking spesific to CNV germline samples
        CNV_GERMLINE_BENCHMARK(
            bench_input.cnv,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(CNV_GERMLINE_BENCHMARK.out.versions)
        ch_reports  = ch_reports.mix(CNV_GERMLINE_BENCHMARK.out.summary_reports)
        ch_summary  = ch_summary.mix(CNV_GERMLINE_BENCHMARK.out.report_multiqc)
    }

    // TODO: SOMATIC BENCHMARKING
    if (params.analysis.contains("somatic")){

        somatic_small = bench_input.snv.mix(bench_input.indel)
        // SOMATIC VARIANT BENCHMARKING
        SMALL_SOMATIC_BENCHMARK(
            somatic_small,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(SMALL_SOMATIC_BENCHMARK.out.versions)
        ch_reports  = ch_reports.mix(SMALL_SOMATIC_BENCHMARK.out.summary_reports)
        ch_summary  = ch_summary.mix(SMALL_SOMATIC_BENCHMARK.out.report_multiqc)

    }

    // compare tool spesfic benchmarks
    COMPARE_BENCHMARK_RESULTS(
        small_evals_ch,
        sv_evals_ch,
        fasta,
        fai
    )
    ch_versions  = ch_versions.mix(COMPARE_BENCHMARK_RESULTS.out.versions)

    // Summarize and plot benchmark statistics
    REPORT_BENCHMARK_STATISTICS(
        ch_reports
    )
    ch_versions = ch_versions.mix(REPORT_BENCHMARK_STATISTICS.out.versions)
    ch_summary  = ch_summary.mix(REPORT_BENCHMARK_STATISTICS.out.report_multiqc)

    // TODO: BENCHMARKING OF CNV
    // https://bioconductor.org/packages/release/bioc/manuals/CNVfilteR/man/CNVfilteR.pdf


    // TODO: TRIO ANALYSIS : MENDELIAN INCONSISTANCE

    // TODO: Compare benchmarking methods!

    //

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(ch_summary)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

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
