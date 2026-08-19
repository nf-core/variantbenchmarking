//
// HAPPY_BENCHMARK: SUBWORKFLOW FOR BENCHMARKING WITH HAPPY
//

include { HAPPY_HAPPY      } from '../../../modules/nf-core/happy/happy/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1    } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2    } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_TRUTH        } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_QUERY        } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_TRUTH_TP } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_TRUTH_FN } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_QUERY_TP } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_QUERY_FP } from '../../../modules/nf-core/bcftools/filter'

workflow HAPPY_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targetsbed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    falsepositive_bed  // reference channel [val(meta), bed]
    stratification_bed // reference channel [val(meta), bed files]
    stratification_tsv // reference channel [val(meta), tsv]

    main:

    tagged_variants = channel.empty()

    input_ch
        .map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
            [ meta, vcf ]
        }
        .set { test_ch }

    input_ch
        .map{ meta, _vcf, _tbi, truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
            [ meta, truth_vcf, _regionsbed, _targets_bed ]
        }
        .set { truth_ch }

    // apply happy method for benchmarking
    HAPPY_HAPPY(
        test_ch.join(truth_ch, failOnDuplicate:true, failOnMismatch:true),
        fasta,
        fai,
        falsepositive_bed,
        stratification_tsv,
        stratification_bed
    )

    // tag meta and collect summary reports
    HAPPY_HAPPY.out.summary_csv
        .map { _meta, csv -> tuple([vartype: params.variant_type] + [benchmark_tool: "happy"], csv) }
        .groupTuple()
        .set{ summary_reports }

    // Subsample TRUTH column from happy results
    BCFTOOLS_VIEW_TRUTH(
        HAPPY_HAPPY.out.vcf.join(HAPPY_HAPPY.out.tbi),
        [],
        [],
        []
        )

    // reheader benchmarking results properly and tag meta
    BCFTOOLS_REHEADER_1(
        BCFTOOLS_VIEW_TRUTH.out.vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )

    // Subsample QUERY column from happy results
    BCFTOOLS_VIEW_QUERY(
        HAPPY_HAPPY.out.vcf.join(HAPPY_HAPPY.out.tbi),
        [],
        [],
        []
    )

    // reheader benchmarking results properly and tag meta
    BCFTOOLS_REHEADER_2(
        BCFTOOLS_VIEW_QUERY.out.vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )

    BCFTOOLS_FILTER_TRUTH_TP(
        BCFTOOLS_REHEADER_1.out.vcf.join(BCFTOOLS_REHEADER_1.out.index)
    )

    BCFTOOLS_FILTER_TRUTH_TP.out.vcf
        .join(BCFTOOLS_FILTER_TRUTH_TP.out.tbi)
        .map { _meta, vcf, index -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "happy"], vcf, index) }
        .set { vcf_tp_comp }

    BCFTOOLS_FILTER_TRUTH_FN(
        BCFTOOLS_REHEADER_1.out.vcf.join(BCFTOOLS_REHEADER_1.out.index)
    )

    BCFTOOLS_FILTER_TRUTH_FN.out.vcf
        .join(BCFTOOLS_FILTER_TRUTH_FN.out.tbi)
        .map { _meta, vcf, index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "happy"], vcf, index) }
        .set { vcf_fn }

    BCFTOOLS_FILTER_QUERY_TP(
        BCFTOOLS_REHEADER_2.out.vcf.join(BCFTOOLS_REHEADER_2.out.index)
    )

    BCFTOOLS_FILTER_QUERY_TP.out.vcf
        .join(BCFTOOLS_FILTER_QUERY_TP.out.tbi)
        .map { _meta, vcf, index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "happy"], vcf, index) }
        .set { vcf_tp_base }

    BCFTOOLS_FILTER_QUERY_FP(
        BCFTOOLS_REHEADER_2.out.vcf.join(BCFTOOLS_REHEADER_2.out.index)
    )

    BCFTOOLS_FILTER_QUERY_FP.out.vcf
        .join(BCFTOOLS_FILTER_QUERY_FP.out.tbi)
        .map { _meta, vcf, index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "happy"], vcf, index) }
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
