/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { RTGTOOLS_BNDEVAL            } from '../modules/nf-core/rtgtools/bndeval'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_variantbenchmarking_pipeline'

//
// SUBWORKFLOWS: Local Subworkflows
//
include { PREPARE_REFERENCES          } from '../subworkflows/local/prepare_references'
include { SUBSAMPLE_VCF_TEST          } from '../subworkflows/local/subsample_vcf_test'
include { PREPARE_VCFS_TRUTH          } from '../subworkflows/local/prepare_vcfs_truth'
include { PREPARE_VCFS_TEST           } from '../subworkflows/local/prepare_vcfs_test'
include { SV_VCF_CONVERSIONS          } from '../subworkflows/local/sv_vcf_conversion'
include { REPORT_VCF_STATISTICS       } from '../subworkflows/local/report_vcf_statistics'
include { CNV_BENCHMARK               } from '../subworkflows/local/cnv_benchmark'
include { SV_BENCHMARK                } from '../subworkflows/local/sv_benchmark'
include { SMALL_BENCHMARK             } from '../subworkflows/local/small_benchmark'
include { REPORT_BENCHMARK_STATISTICS } from '../subworkflows/local/report_benchmark_statistics'
include { COMPARE_BENCHMARK_RESULTS   } from '../subworkflows/local/compare_benchmark_results'
include { INTERSECT_STATISTICS        } from '../subworkflows/local/intersect_statistics'
include { CONCORDANCE_ANALYSIS        } from '../subworkflows/local/concordance_analysis'
include { ENSEMBLE_TEST_VCFS          } from '../subworkflows/local/ensemble_test_vcfs'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTBENCHMARKING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // To gather all QC reports for Multiqc
    ch_versions      = channel.empty()
    ch_multiqc_files = channel.empty()

    // To gather benchmark reports
    ch_reports       = channel.empty()

    //// create reference channels ////

    fasta       = channel.fromPath(params.fasta, checkIfExists: true)
                    .map{ fasta -> tuple([id: fasta.getSimpleName()], fasta) }.collect()
    fai         = channel.fromPath(params.fai, checkIfExists: true)
                    .map{ fai -> tuple([id: fai.getSimpleName()], fai) }.collect()

    //// check Truth Files ////
    truth_ch        = params.truth_vcf ? channel.fromPath(params.truth_vcf, checkIfExists: true)
                                            .map{ vcf -> tuple([id: params.truth_id, vartype:params.variant_type], vcf) }.collect()
                                        : channel.empty()

    regions_bed_ch = params.regions_bed ? channel.fromPath(params.regions_bed, checkIfExists: true).collect()
                                        : channel.empty()
    targets_bed_ch = params.targets_bed ? channel.fromPath(params.targets_bed, checkIfExists: true).collect()
                                        : channel.empty()



    if (params.method ==~ /.*(?:truvari|svanalyzer|happy|sompy|rtgtools|wittyer).*/) {
        // Note: concordance analysis does not require truth files
        if (params.ensemble_truth){
            log.warn "[nf-core/variantbenchmarking] WARN: --truth_id will be treated as 'truth', meaning, ensembled truth file will be named as 'truth'"
            } else
            {
                if (!params.truth_vcf || !params.truth_id){
                error "[nf-core/variantbenchmarking] ERROR:Please specify --truth_id and --truth_vcf to perform benchmarking analysis"
            }
        }
    }
    if (params.method ==~ /.*(?:intersect).*/) {
        if (!params.regions_bed || !params.truth_vcf){
            error "[nf-core/variantbenchmarking] ERROR: Please specify --regions_bed to perform intersect analysis"
        }
    }

    if (params.preprocess?.contains("filter_contigs") && !params.genome) {
        error "[nf-core/variantbenchmarking] ERROR:Please specify --genome when using filter_contigs for --preprocessing"
    }

    // Optional files for Happy or Sompy
    ambiguous_beds      = params.ambiguous_beds     ? channel.fromPath(params.ambiguous_beds, checkIfExists: true).map{ bed -> tuple([id: "ambiguous"], bed) }.collect()
                                                    : channel.of([[id: "ambiguous"],[]]).collect()
    if (params.stratification_bed && params.stratification_tsv){
        stratification_bed  = channel.fromPath(params.stratification_bed, checkIfExists: true, type: 'dir').map{ bed -> tuple([id: "stratification"], bed) }.collect()
        stratification_tsv  = channel.fromPath(params.stratification_tsv, checkIfExists: true).map{ tsv -> tuple([id: "stratification"], tsv) }.collect()
    }else{
        stratification_bed  = channel.of([[id: "stratification"],[]]).collect()
        stratification_tsv  = channel.of([[id: "stratification"],[]]).collect()
    }

    // SDF file for RTG-tools eval
    sdf             = params.sdf        ? channel.fromPath(params.sdf, checkIfExists: true).map{ sdf -> tuple([id: sdf.getSimpleName()], sdf) }.collect()
                                        : channel.empty()

    if (params.rename_chr){
        rename_chr = channel.fromPath(params.rename_chr, checkIfExists: true).map{ txt -> tuple([id: txt.getSimpleName()], txt) }.collect()
        if (!params.genome){
            error "[nf-core/variantbenchmarking] ERROR:Please specify --genome to fix chromosome prefix"
        }

    }else{
        if (params.genome == "GRCh38"){
            rename_chr = channel.fromPath("${projectDir}/assets/rename_contigs/grch37_grch38.txt", checkIfExists: true).map{ txt -> tuple([id: txt.getSimpleName()], txt) }.collect()
        }
        else if(params.genome == "GRCh37")
        {
            rename_chr = channel.fromPath("${projectDir}/assets/rename_contigs/grch38_grch37.txt", checkIfExists: true).map{ txt -> tuple([id: txt.getSimpleName()], txt) }.collect()
        }
        else{
            rename_chr = channel.empty()
        }
    }

    // read chain file, liftover genome and rename chr files if liftover is true
    if (params.liftover){

        if (params.chain){
            chain       = channel.fromPath(params.chain, checkIfExists: true).map{ bed -> tuple([id: bed.getSimpleName()], bed) }.collect()
        }else{
            error "[nf-core/variantbenchmarking] ERROR:Please specify --chain to process liftover of the files"
        }
        // if dictionary file is missing PICARD_CREATESEQUENCEDICTIONARY will create one
        dictionary      = params.dictionary ? channel.fromPath(params.dictionary, checkIfExists: true).map{ dict -> tuple([id: dict.getSimpleName()], dict) }.collect()
                                            : channel.empty()
    }else{
        chain           = channel.empty()
        dictionary      = channel.empty()
    }

    //prepare references and libraries
    PREPARE_REFERENCES(
        fasta,
        dictionary,
        sdf
    )
    dictionary  = PREPARE_REFERENCES.out.dictionary
    sdf         = PREPARE_REFERENCES.out.sdf

    // PREPROCESSES

    // subsample multisample vcf if necessary, filter out cases without test vcf (only regions)
    ch_samplesheet.branch{ input ->
            def meta = input[0]
            def vcf = input[1]
            multisample: meta.subsample
            singlesample : vcf
            other: false}
        .set{sample}

    out_vcf_ch  = channel.empty()

    SUBSAMPLE_VCF_TEST(
        sample.multisample.map{meta, vcf, _bed -> [meta, vcf]}
    )
    vcf_ch  = out_vcf_ch.mix(SUBSAMPLE_VCF_TEST.out.vcf_ch,
                                sample.singlesample.map{meta, vcf, _bed -> [meta, vcf]})

    if (params.variant_type == "structural" || params.variant_type == "copynumber" ){
        // Standardize SV VCFs, tool specific modifications
        SV_VCF_CONVERSIONS(
            vcf_ch,
            fai
        )
        vcf_ch      = SV_VCF_CONVERSIONS.out.vcf_ch.map{it -> tuple(it[0], it[1])}
    }
    // Prepare and normalize input vcfs
    PREPARE_VCFS_TEST(
        vcf_ch,
        fasta,
        fai,
        chain,
        rename_chr,
        dictionary
    )

    // Ensemble and prepare truth file using input VCFs if ensemble_truth approach is choosen
    if (params.ensemble_truth){
        ENSEMBLE_TEST_VCFS(
            PREPARE_VCFS_TEST.out.vcf_ch,
            fasta,
            fai
        )
        truth_ch    = ENSEMBLE_TEST_VCFS.out.truth_vcf
    } else {
        // Prepare and normalize truth vcf provided
        PREPARE_VCFS_TRUTH(
            truth_ch,
            regions_bed_ch,
            targets_bed_ch,
            fasta,
            fai,
            chain,
            rename_chr,
            dictionary
        )
        regions_bed_ch = PREPARE_VCFS_TRUTH.out.regions_bed_ch
        targets_bed_ch = PREPARE_VCFS_TRUTH.out.targets_bed_ch
        truth_ch       = PREPARE_VCFS_TRUTH.out.vcf_ch
    }

    // VCF REPORTS AND STATS

    // get statistics for normalized input files
    REPORT_VCF_STATISTICS(
        PREPARE_VCFS_TEST.out.vcf_ch.mix(truth_ch)
    )
    ch_multiqc_files = ch_multiqc_files.mix(REPORT_VCF_STATISTICS.out.ch_stats)


    // If intersect is in the methods, perform bedtools intersect to region files given
    ch_samplesheet.branch{
        def _meta = it[0]
        def regions_file = it[2]
        regions : regions_file
        other: false}
        .set{intersect}

    if (params.method.contains("intersect")){

        INTERSECT_STATISTICS(
            intersect.regions.mix(PREPARE_VCFS_TEST.out.vcf_ch),
            regions_bed_ch
        )
        ch_reports = ch_reports.mix(INTERSECT_STATISTICS.out.summary_reports)
    }
    evals_ch     = channel.empty()
    evals_csv_ch = channel.empty()

    // Concordance analysis can only be performed small variants for now
    if (params.method.contains("concordance") && (params.variant_type ==~ /.*(?:small|snv|indel).*/)){
      CONCORDANCE_ANALYSIS(
            PREPARE_VCFS_TEST.out.vcf_ch,
            regions_bed_ch,
            fasta,
            fai,
            dictionary
        )
        ch_reports       = ch_reports.mix(CONCORDANCE_ANALYSIS.out.summary_reports)
        evals_ch         = evals_ch.mix(CONCORDANCE_ANALYSIS.out.tagged_variants)

    }else if (params.method.contains("concordance")){
        error "[nf-core/variantbenchmarking] ERROR:Concordance analysis can only be performed small (snv and indels) variants"
    }

    // Prepare benchmark channel
    PREPARE_VCFS_TEST.out.vcf_ch.combine(truth_ch)
        .combine(regions_bed_ch.ifEmpty([[]]))
        .combine(targets_bed_ch.ifEmpty([[]]))
        .map{ test_meta, test_vcf, test_tbi, _truth_meta, truth_vcf, truth_tbi, regions_bed, targets_bed  ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, regions_bed, targets_bed ]}
        .set{bench}


    if (params.variant_type == "structural"){
        SV_BENCHMARK(
            bench,
            fasta,
            fai
        )
        ch_reports       = ch_reports.mix(SV_BENCHMARK.out.summary_reports)
        evals_ch         = evals_ch.mix(SV_BENCHMARK.out.tagged_variants)
        ch_multiqc_files = ch_multiqc_files.mix(SV_BENCHMARK.out.logs.map{ _meta, log -> log })

        if (params.method.contains("rtgtools")){
            // running bndeval might require svdecompose
            RTGTOOLS_BNDEVAL(
                bench.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
                    [ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed ]
                }
            )

            ch_reports  = ch_reports.mix(RTGTOOLS_BNDEVAL.out.summary
                                                    .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], file) }
                                                    .groupTuple())

            evals_ch = evals_ch.mix(
                RTGTOOLS_BNDEVAL.out.fn_vcf.map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "rtgtools"], file) },
                RTGTOOLS_BNDEVAL.out.fp_vcf.map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "rtgtools"], file) },
                RTGTOOLS_BNDEVAL.out.baseline_vcf.map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "rtgtools"], file) },
                RTGTOOLS_BNDEVAL.out.tp_vcf.map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "rtgtools"], file) }
            )
        }
    }
    if ( params.variant_type == "copynumber"){
        CNV_BENCHMARK(
            bench,
            fasta,
            fai
        )
        ch_reports       = ch_reports.mix(CNV_BENCHMARK.out.summary_reports)
        evals_ch         = evals_ch.mix(CNV_BENCHMARK.out.tagged_variants)
        ch_multiqc_files = ch_multiqc_files.mix(CNV_BENCHMARK.out.logs.map{ _meta, log -> log })
    }

    if (params.variant_type == "small" || params.variant_type == "snv" || params.variant_type == "indel"){
        // Benchmarking small variants
        SMALL_BENCHMARK(
            bench,
            fasta,
            fai,
            sdf,
            stratification_bed,
            stratification_tsv,
            ambiguous_beds
        )
        ch_reports       = ch_reports.mix(SMALL_BENCHMARK.out.summary_reports)
        evals_ch         = evals_ch.mix(SMALL_BENCHMARK.out.tagged_variants)
        evals_csv_ch     = evals_csv_ch.mix(SMALL_BENCHMARK.out.tagged_variants_csv)
    }


    // Compare tagged variants
    COMPARE_BENCHMARK_RESULTS(
        evals_ch,
        evals_csv_ch,
        fasta,
        fai
    )
    ch_multiqc_files = ch_multiqc_files.mix(COMPARE_BENCHMARK_RESULTS.out.ch_plots.flatten())

    // Summarize and plot benchmark statistics
    REPORT_BENCHMARK_STATISTICS(
        ch_reports,
        evals_ch,
        evals_csv_ch
    )
    ch_multiqc_files = ch_multiqc_files.mix(REPORT_BENCHMARK_STATISTICS.out.ch_plots.flatten())
    ch_multiqc_files = ch_multiqc_files.mix(REPORT_BENCHMARK_STATISTICS.out.merged_reports.map{ _meta, report -> report })
    ch_multiqc_files = ch_multiqc_files.mix(ch_reports.map{ _meta, report -> report }.flatten())

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'variantbenchmarking_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))


    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'variantbenchmarking'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
