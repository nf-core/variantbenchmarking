//
// RTGTOOLS_BENCHMARK: SUBWORKFLOW FOR RTGTOOLS_BENCHMARKING
//

include { RTGTOOLS_FORMAT  } from '../../../modules/nf-core/rtgtools/format/main'
include { RTGTOOLS_VCFEVAL } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1    } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2    } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_3    } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_4    } from '../../../modules/local/bcftools/reheader'

workflow RTGTOOLS_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targetsbed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    sdf                // reference channel [val(meta), sdf]

    main:

    versions        = Channel.empty()
    tagged_variants = Channel.empty()

    if (!params.sdf){

        // Use rtgtools format to generate sdf file if necessary
        RTGTOOLS_FORMAT(
            fasta.map { meta, file -> [ meta, file, [], [] ] }
        )
        versions = versions.mix(RTGTOOLS_FORMAT.out.versions)
        sdf = RTGTOOLS_FORMAT.out.sdf
    }

    // apply rtgtools eval method
    RTGTOOLS_VCFEVAL(
        input_ch,
        sdf
    )
    versions = versions.mix(RTGTOOLS_VCFEVAL.out.versions.first())

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
    versions = versions.mix(BCFTOOLS_REHEADER_1.out.versions.first())

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
    versions = versions.mix(BCFTOOLS_REHEADER_2.out.versions)

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
    versions = versions.mix(BCFTOOLS_REHEADER_3.out.versions)

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

    emit:
    summary_reports // channel: [val(meta), reports]
    tagged_variants // channel: [val(meta), vcfs]
    versions        // channel: [versions.yml]

}
