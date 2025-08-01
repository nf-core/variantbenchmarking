/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_variantbenchmarking_pipeline'

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
include { INTERSECT_STATISTICS        } from '../subworkflows/local/intersect_statistics'
include { BND_BENCHMARK               } from '../subworkflows/local/bnd_benchmark'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTBENCHMARKING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    // To gather all QC reports for Multiqc
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // To gather benchmark reports
    ch_reports       = Channel.empty()

    //// create reference channels ////

    fasta       = Channel.fromPath(params.fasta, checkIfExists: true)
                    .map{ fasta -> tuple([id: fasta.getSimpleName()], fasta) }.collect()
    fai         = Channel.fromPath(params.fai, checkIfExists: true)
                    .map{ fai -> tuple([id: fai.getSimpleName()], fai) }.collect()

    //// check Truth Files ////

    if (params.truth_vcf || params.regions_bed){
        truth_ch        = params.truth_vcf ? Channel.fromPath(params.truth_vcf, checkIfExists: true)
                                                    .map{ vcf -> tuple([id: params.truth_id, vartype:params.variant_type], vcf) }.collect()
                                                    : Channel.empty()

        regions_bed_ch = params.regions_bed ? Channel.fromPath(params.regions_bed, checkIfExists: true).collect()
                                                    : Channel.empty()
        targets_bed_ch = params.targets_bed ? Channel.fromPath(params.targets_bed, checkIfExists: true).collect()
                                                    : Channel.empty()
    }else{
        log.error "Please specify params.truth_id and params.truth_vcf or params.regions_bed to perform benchmarking analysis"
        exit 1
    }

    // Optional files for Happy or Sompy
    falsepositive_bed   = params.falsepositive_bed  ? Channel.fromPath(params.falsepositive_bed, checkIfExists: true).map{ bed -> tuple([id: "falsepositive"], bed) }.collect()
                                                    : Channel.of([[id: "falsepositive"],[]]).collect()
    ambiguous_beds      = params.ambiguous_beds     ? Channel.fromPath(params.ambiguous_beds, checkIfExists: true).map{ bed -> tuple([id: "ambiguous"], bed) }.collect()
                                                    : Channel.of([[id: "ambiguous"],[]]).collect()
    if (params.stratification_bed && params.stratification_tsv){
        stratification_bed  = Channel.fromPath(params.stratification_bed, checkIfExists: true, type: 'dir').map{ bed -> tuple([id: "stratification"], bed) }.collect()
        stratification_tsv  = Channel.fromPath(params.stratification_tsv, checkIfExists: true).map{ tsv -> tuple([id: "stratification"], tsv) }.collect()
    }else{
        stratification_bed  = Channel.of([[id: "stratification"],[]]).collect()
        stratification_tsv  = Channel.of([[id: "stratification"],[]]).collect()
    }

    // SDF file for RTG-tools eval
    sdf             = params.sdf        ? Channel.fromPath(params.sdf, checkIfExists: true).map{ sdf -> tuple([id: sdf.getSimpleName()], sdf) }.collect()
                                        : Channel.empty()

    if (params.rename_chr){
        rename_chr = Channel.fromPath(params.rename_chr, checkIfExists: true).map{ txt -> tuple([id: txt.getSimpleName()], txt) }.collect()
        if (!params.genome){
            log.error "Please specify params.genome to fix chromosome prefix"
            exit 1
        }

    }else{
        if (params.genome == "GRCh38"){
            rename_chr = Channel.fromPath("${projectDir}/assets/rename_contigs/grch37_grch38.txt", checkIfExists: true).map{ txt -> tuple([id: txt.getSimpleName()], txt) }.collect()
        }
        else if(params.genome == "GRCh37")
        {
            rename_chr = Channel.fromPath("${projectDir}/assets/rename_contigs/grch38_grch37.txt", checkIfExists: true).map{ txt -> tuple([id: txt.getSimpleName()], txt) }.collect()
        }
        else{
            rename_chr = Channel.empty()
        }
    }

    // read chain file, liftover genome and rename chr files if liftover is true
    if (params.liftover){

        if (params.chain){
            chain           = Channel.fromPath(params.chain, checkIfExists: true).map{ bed -> tuple([id: bed.getSimpleName()], bed) }.collect()
        }else{
            log.error "Please specify params.chain to process liftover of the files"
            exit 1
        }
        // if dictionary file is missing PICARD_CREATESEQUENCEDICTIONARY will create one
        dictionary      = params.dictionary ? Channel.fromPath(params.dictionary, checkIfExists: true).map{ dict -> tuple([id: dict.getSimpleName()], dict) }.collect()                                           : Channel.empty()
    }else{
        chain           = Channel.empty()
        dictionary      = Channel.empty()
    }

    if(params.method.contains("intersect")){
        if(!params.regions_bed){
            log.error "Regions BED is required for intersection analysis"
            exit 1
        }
    }

    // PREPROCESSES

    // subsample multisample vcf if necessary, filter out cases without test vcf (only regions)

    ch_samplesheet.branch{
            def meta = it[0]
            def vcf = it[1]
            multisample: meta.subsample
            singlesample : vcf
            other: false}
        .set{sample}

    out_vcf_ch  = Channel.empty()

    SUBSAMPLE_VCF_TEST(
        sample.multisample.map{meta, vcf, bed -> [meta, vcf]}
    )
    ch_versions = ch_versions.mix(SUBSAMPLE_VCF_TEST.out.versions)
    vcf_ch  = out_vcf_ch.mix(SUBSAMPLE_VCF_TEST.out.vcf_ch,
                                sample.singlesample.map{meta, vcf, bed -> [meta, vcf]})

    if (params.variant_type == "structural" ){
        // Standardize SV VCFs, tool specific modifications
        SV_VCF_CONVERSIONS(
            vcf_ch,
            fasta,
            fai
        )
        ch_versions = ch_versions.mix(SV_VCF_CONVERSIONS.out.versions)
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
    ch_versions = ch_versions.mix(PREPARE_VCFS_TEST.out.versions)

    // Prepare and normalize truth vcfs
    PREPARE_VCFS_TRUTH(
        truth_ch,
        regions_bed_ch,
        fasta,
        fai,
        chain,
        rename_chr,
        dictionary
    )
    regions_bed_ch = PREPARE_VCFS_TRUTH.out.high_conf_ch
    ch_versions    = ch_versions.mix(PREPARE_VCFS_TRUTH.out.versions)

    // VCF REPORTS AND STATS

    // get statistics for normalized input files
    REPORT_VCF_STATISTICS(
        PREPARE_VCFS_TEST.out.vcf_ch.mix(PREPARE_VCFS_TRUTH.out.vcf_ch)
    )
    ch_versions       = ch_versions.mix(REPORT_VCF_STATISTICS.out.versions)
    ch_multiqc_files  = ch_multiqc_files.mix(REPORT_VCF_STATISTICS.out.ch_stats)


    // If intersect is in the methods, perform bedtools intersect to region files given
    ch_samplesheet.branch{
        def meta = it[0]
        def regions_file = it[2]
        regions : regions_file
        other: false}
        .set{intersect}

    if (params.method.contains("intersect")){

        INTERSECT_STATISTICS(
            intersect.regions.mix(PREPARE_VCFS_TEST.out.vcf_ch),
            regions_bed_ch
        )
        ch_versions      = ch_versions.mix(INTERSECT_STATISTICS.out.versions)
        ch_reports       = ch_reports.mix(INTERSECT_STATISTICS.out.summary_reports)

    }

    // Prepare benchmark channel
    PREPARE_VCFS_TEST.out.vcf_ch.combine(PREPARE_VCFS_TRUTH.out.vcf_ch)
        .combine(regions_bed_ch.ifEmpty([[]]))
        .combine(targets_bed_ch.ifEmpty([[]]))
        .map{ test_meta, test_vcf, test_tbi, _truth_meta, truth_vcf, truth_tbi, regions_bed, targets_bed  ->
                    [ test_meta, test_vcf, test_tbi, truth_vcf, truth_tbi, regions_bed, targets_bed ]}
        .set{bench}

    evals_ch     = Channel.empty()
    evals_csv_ch = Channel.empty()


    if (params.variant_type == "structural" || params.variant_type == "copynumber"){
        // Perform SV benchmarking - for now it also works for somatic cases!
        // this part will be changed!
        SV_GERMLINE_BENCHMARK(
            bench,
            fasta,
            fai
        )
        ch_versions      = ch_versions.mix(SV_GERMLINE_BENCHMARK.out.versions)
        ch_reports       = ch_reports.mix(SV_GERMLINE_BENCHMARK.out.summary_reports)
        evals_ch         = evals_ch.mix(SV_GERMLINE_BENCHMARK.out.tagged_variants)

        if (params.method.contains("bndeval")){
            // running bndeval might require svdecompose
            BND_BENCHMARK(
                bench,
                fai
            )

            ch_versions      = ch_versions.mix(BND_BENCHMARK.out.versions)
            ch_reports       = ch_reports.mix(BND_BENCHMARK.out.summary_reports)
            evals_ch         = evals_ch.mix(BND_BENCHMARK.out.tagged_variants)
        }
    }

    if (params.analysis.contains("germline")){

        if (params.variant_type == "small" | params.variant_type == "snv" | params.variant_type == "indel"){
            // Benchmarking specific to small germline samples
            SMALL_GERMLINE_BENCHMARK(
                bench,
                fasta,
                fai,
                sdf,
                falsepositive_bed,
                stratification_bed,
                stratification_tsv
            )
            ch_versions      = ch_versions.mix(SMALL_GERMLINE_BENCHMARK.out.versions)
            ch_reports       = ch_reports.mix(SMALL_GERMLINE_BENCHMARK.out.summary_reports)
            evals_ch         = evals_ch.mix(SMALL_GERMLINE_BENCHMARK.out.tagged_variants)
        }
    }

    // SOMATIC BENCHMARKING
    if (params.analysis.contains("somatic")){

        if (params.variant_type == "snv" | params.variant_type == "indel"){
            // SOMATIC VARIANT BENCHMARKING
            SMALL_SOMATIC_BENCHMARK(
                bench,
                fasta,
                fai,
                sdf,
                falsepositive_bed,
                ambiguous_beds
            )
            ch_versions      = ch_versions.mix(SMALL_SOMATIC_BENCHMARK.out.versions)
            ch_reports       = ch_reports.mix(SMALL_SOMATIC_BENCHMARK.out.summary_reports)
            evals_ch         = evals_ch.mix(SMALL_SOMATIC_BENCHMARK.out.tagged_variants)
            evals_csv_ch     = evals_csv_ch.mix(SMALL_SOMATIC_BENCHMARK.out.tagged_variants_csv)
        }

    }

    // compare tool spesfic benchmarks
    COMPARE_BENCHMARK_RESULTS(
        evals_ch,
        evals_csv_ch,
        fasta,
        fai
    )
    ch_versions  = ch_versions.mix(COMPARE_BENCHMARK_RESULTS.out.versions)

    // Summarize and plot benchmark statistics
    REPORT_BENCHMARK_STATISTICS(
        ch_reports
    )
    ch_versions      = ch_versions.mix(REPORT_BENCHMARK_STATISTICS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'variantbenchmarking_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) :Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml',sort: true))
    ch_multiqc_files                      = ch_multiqc_files.mix(REPORT_BENCHMARK_STATISTICS.out.ch_plots)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions            = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
