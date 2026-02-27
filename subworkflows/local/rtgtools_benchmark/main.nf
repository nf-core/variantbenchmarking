//
// RTGTOOLS_BENCHMARK: SUBWORKFLOW FOR RTGTOOLS_BENCHMARKING
//

include { RTGTOOLS_VCFEVAL } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1    } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2    } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_3    } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_4    } from '../../../modules/nf-core/bcftools/reheader'

workflow RTGTOOLS_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targetsbed ]
    fai                // reference channel [val(meta), ref.fa.fai]
    sdf                // reference channel [val(meta), sdf]

    main:

    tagged_variants = channel.empty()

    // apply rtgtools eval method
    RTGTOOLS_VCFEVAL(
        input_ch,
        sdf
    )

    // collect summary reports
    RTGTOOLS_VCFEVAL.out.summary
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], file) }
        .groupTuple()
        .set{ summary_reports }

    // reheader benchmarking results properly and tag meta
    BCFTOOLS_REHEADER_1(
        RTGTOOLS_VCFEVAL.out.fn_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )

    BCFTOOLS_REHEADER_1.out.vcf
        .join(BCFTOOLS_REHEADER_1.out.index)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "rtgtools"], file, index) }
        .set { vcf_fn }

    BCFTOOLS_REHEADER_2(
        RTGTOOLS_VCFEVAL.out.fp_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )

    BCFTOOLS_REHEADER_2.out.vcf
        .join(BCFTOOLS_REHEADER_2.out.index)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "rtgtools"], file, index) }
        .set { vcf_fp }

    BCFTOOLS_REHEADER_3(
        RTGTOOLS_VCFEVAL.out.baseline_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )

    BCFTOOLS_REHEADER_3.out.vcf
        .join(BCFTOOLS_REHEADER_3.out.index)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "rtgtools"], file, index) }
        .set { vcf_tp_base }

    BCFTOOLS_REHEADER_4(
        RTGTOOLS_VCFEVAL.out.tp_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )

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

    emit:
    summary_reports // channel: [val(meta), reports]
    tagged_variants // channel: [val(meta), vcfs]

}
