//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

include { HAPPY_SOMPY          } from '../../../modules/nf-core/happy/sompy'
include { RTGTOOLS_FORMAT      } from '../../../modules/nf-core/rtgtools/format'
include { SPLIT_SOMPY_FEATURES } from '../../../modules/local/custom/split_sompy_features'
include { RTGTOOLS_VCFEVAL  as RTGTOOLS_VCFEVAL_SOMATIC  } from '../../../modules/nf-core/rtgtools/vcfeval'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1       } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2       } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_3       } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_4       } from '../../../modules/local/bcftools/reheader'

workflow SMALL_SOMATIC_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index,  regionsbed, targets_bed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    sdf                // reference channel [val(meta), sdf]
    falsepositive_bed  // reference channel [val(meta), bed]
    ambiguous_beds     // reference channel [val(meta), bed]

    main:

    versions            = Channel.empty()
    summary_reports     = Channel.empty()
    tagged_variants     = Channel.empty()
    tagged_variants_csv = Channel.empty()

    if (params.method.contains('sompy')){
        // apply sompy for small somatic variant benchmarking
        HAPPY_SOMPY(
            input_ch.map{meta, test, test_index, truth, truth_index, regions, target -> [meta, test, truth, regions, target]},
            fasta,
            fai,
            falsepositive_bed,
            ambiguous_beds,
            [[],[]]
        )
        versions = versions.mix(HAPPY_SOMPY.out.versions.first())

        HAPPY_SOMPY.out.stats
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "sompy"], file) }
            .groupTuple()
            .set{ report }
        summary_reports = summary_reports.mix(report)

        SPLIT_SOMPY_FEATURES(
            HAPPY_SOMPY.out.features
        )
        versions = versions.mix(SPLIT_SOMPY_FEATURES.out.versions)

        SPLIT_SOMPY_FEATURES.out.TP
            .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "TP"] + [id: "sompy"], file) }
            .set{tp_vars}

        SPLIT_SOMPY_FEATURES.out.FP
            .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "sompy"], file) }
            .set{fp_vars}

        SPLIT_SOMPY_FEATURES.out.FN
            .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "sompy"], file) }
            .set{fn_vars}

        tagged_variants_csv = tagged_variants_csv
                                .mix(tp_vars)
                                .mix(fp_vars)
                                .mix(fn_vars)
    }

    if (params.method.contains('rtgtools')){

        if (!params.sdf){

            // Use rtgtools format to generate sdf file if necessary
            RTGTOOLS_FORMAT(
                fasta.map { meta, file -> [ meta, file, [], [] ] }
            )
            versions = versions.mix(RTGTOOLS_FORMAT.out.versions)
            sdf = RTGTOOLS_FORMAT.out.sdf

        }

        // apply rtgtools eval method
        RTGTOOLS_VCFEVAL_SOMATIC(
            input_ch,
            sdf
        )
        versions = versions.mix(RTGTOOLS_VCFEVAL_SOMATIC.out.versions.first())

        // collect summary reports
        RTGTOOLS_VCFEVAL_SOMATIC.out.summary
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], file) }
            .groupTuple()
            .set{ report }

        summary_reports = summary_reports.mix(report)

        // reheader benchmarking results properly and tag meta
        BCFTOOLS_REHEADER_1(
            RTGTOOLS_VCFEVAL_SOMATIC.out.fn_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_1.out.versions.first())

        BCFTOOLS_REHEADER_1.out.vcf
            .join(BCFTOOLS_REHEADER_1.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "rtgtools"], file, index) }
            .set { vcf_fn }

        BCFTOOLS_REHEADER_2(
            RTGTOOLS_VCFEVAL_SOMATIC.out.fp_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_2.out.versions)

        BCFTOOLS_REHEADER_2.out.vcf
            .join(BCFTOOLS_REHEADER_2.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "rtgtools"], file, index) }
            .set { vcf_fp }

        BCFTOOLS_REHEADER_3(
            RTGTOOLS_VCFEVAL_SOMATIC.out.baseline_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_3.out.versions)

        BCFTOOLS_REHEADER_3.out.vcf
            .join(BCFTOOLS_REHEADER_3.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "rtgtools"], file, index) }
            .set { vcf_tp_base }

        BCFTOOLS_REHEADER_4(
            RTGTOOLS_VCFEVAL_SOMATIC.out.tp_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_4.out.versions)

        BCFTOOLS_REHEADER_4.out.vcf
            .join(BCFTOOLS_REHEADER_4.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "rtgtools"], file, index) }
            .set { vcf_tp_comp }

        tagged_variants = tagged_variants.mix(
            vcf_fn,
            vcf_fp,
            vcf_tp_base,
            vcf_tp_comp
        )
    }

    emit:
    summary_reports     // channel: [val(meta), reports]
    tagged_variants     // channel: [val(meta), vcfs]
    tagged_variants_csv // channel: [val(meta), csvs]
    versions            // channel: [versions.yml]
}
