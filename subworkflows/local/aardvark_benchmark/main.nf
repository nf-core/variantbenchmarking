//
// AARDVARK_BENCHMARK: SUBWORKFLOW FOR BENCHMARKING WITH AARDVARK
//

include { AARDVARK_COMPARE  } from '../../../modules/nf-core/aardvark/compare'
include { BCFTOOLS_REHEADER } from '../../../modules/nf-core/bcftools/reheader'
include { TABIX_TABIX       } from '../../../modules/nf-core/tabix/tabix'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_TRUTH_TP } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_TRUTH_FN } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_QUERY_TP } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_QUERY_FP } from '../../../modules/nf-core/bcftools/filter'

workflow AARDVARK_BENCHMARK {

    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targetsbed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    stratification_bed // reference channel [val(meta), bed files]
    stratification_tsv // reference channel [val(meta), tsv]

    main:

    tagged_variants = channel.empty()

    AARDVARK_COMPARE(
        input_ch,
        fasta,
        stratification_bed,
        stratification_tsv
    )

    AARDVARK_COMPARE.out.summary
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "aardvark"], file) }
            .groupTuple()
            .map { meta, files -> tuple(meta, files.flatten()) }
            .set { summary_reports }

    // Filter TP/FN from labelled_truth
    // reheader truth vcf with query names to enable comparisons better for plotting
    BCFTOOLS_REHEADER(
         AARDVARK_COMPARE.out.labelled_truth.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )

    BCFTOOLS_FILTER_TRUTH_TP(
        BCFTOOLS_REHEADER.out.vcf.join(BCFTOOLS_REHEADER.out.index)
    )

    BCFTOOLS_FILTER_TRUTH_TP.out.vcf
        .join(BCFTOOLS_FILTER_TRUTH_TP.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "aardvark"], file, index) }
        .set { vcf_tp_comp }

    BCFTOOLS_FILTER_TRUTH_FN(
        BCFTOOLS_REHEADER.out.vcf.join(BCFTOOLS_REHEADER.out.index)
    )

    BCFTOOLS_FILTER_TRUTH_FN.out.vcf
        .join(BCFTOOLS_FILTER_TRUTH_FN.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "aardvark"], file, index) }
        .set { vcf_fn }

    // Filter TP/FP from labelled_query
    TABIX_TABIX(
        AARDVARK_COMPARE.out.labelled_query
    )

    BCFTOOLS_FILTER_QUERY_TP(
        AARDVARK_COMPARE.out.labelled_query.join(TABIX_TABIX.out.index)
    )

    BCFTOOLS_FILTER_QUERY_TP.out.vcf
        .join(BCFTOOLS_FILTER_QUERY_TP.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "aardvark"], file, index) }
        .set { vcf_tp_base }

    BCFTOOLS_FILTER_QUERY_FP(
        AARDVARK_COMPARE.out.labelled_query.join(TABIX_TABIX.out.index)
    )

    BCFTOOLS_FILTER_QUERY_FP.out.vcf
        .join(BCFTOOLS_FILTER_QUERY_FP.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "aardvark"], file, index) }
        .set { vcf_fp }

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
