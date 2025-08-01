//
// INTERSECTION ANALYSIS OF BED FILES
//

include { BEDTOOLS_INTERSECT_BENCH } from '../../../modules/local/custom/bedtools_intersect_bench'
include { SVTK_VCF2BED             } from '../../../modules/nf-core/svtk/vcf2bed'
include { BEDOPS_CONVERT2BED       } from '../../../modules/nf-core/bedops/convert2bed'
include { TABIX_BGZIP              } from '../../../modules/nf-core/tabix/bgzip'

workflow INTERSECT_STATISTICS {
    take:
    test             // channel: [val(meta), vcf, regions]
    truth_regions    // channel: [truth bed]

    main:

    versions        = Channel.empty()

    test.branch{
            def vcf_file = it[1]
            def regions_file = it[2]
            regions : regions_file.extension != "tbi"
            vcf : vcf_file
            other: false}
            .set{test_samples}


    test_beds_ch = test_samples.regions.map{meta, vcf, bed -> [meta, bed]}

    // convert VCF files to BED format
    if (params.variant_type == "structural" || params.variant_type == "copynumber"){
        SVTK_VCF2BED(
            test_samples.vcf
        )

        // collect summary reports
        SVTK_VCF2BED.out.bed
            .map { meta, file -> tuple(meta + [converted: true], file) }
            .set{ converted_beds }
        test_beds_ch = test_beds_ch.mix(converted_beds)
        versions  = versions.mix(SVTK_VCF2BED.out.versions)

    }

    if (params.variant_type == "small" || params.variant_type == "snv" || params.variant_type == "indel"){

        // unzip vcf.gz file
        TABIX_BGZIP(
            test_samples.vcf.map { meta, vcf, _tbi->[ meta, vcf ]}
        )
        versions  = versions.mix(TABIX_BGZIP.out.versions)

        BEDOPS_CONVERT2BED(
            TABIX_BGZIP.out.output
        )
        // collect summary reports
        BEDOPS_CONVERT2BED.out.bed
            .map { meta, file -> tuple(meta + [converted: true], file) }
            .set{ converted_beds }
        test_beds_ch = test_beds_ch.mix(converted_beds)
        versions  = versions.mix(BEDOPS_CONVERT2BED.out.versions)
    }

    test_beds_ch
            .combine(truth_regions)
            .map{test_meta, testbed, truthbed -> [test_meta, truthbed, testbed]}
            .set{intersect_ch}

    // Intersect bed files and gather statistics
    BEDTOOLS_INTERSECT_BENCH(
        intersect_ch
    )
    versions      = versions.mix(BEDTOOLS_INTERSECT_BENCH.out.versions)

    // collect summary reports
    BEDTOOLS_INTERSECT_BENCH.out.summary
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "intersect"], file) }
        .groupTuple()
        .set{ summary_reports }

    emit:
    versions         // channel: [versions.yml]
    summary_reports  // channel: [meta, summary_report.csv]
}
