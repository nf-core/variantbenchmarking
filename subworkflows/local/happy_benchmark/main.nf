//
// HAPPY_BENCHMARK: SUBWORKFLOW FOR BENCHMARKING WITH HAPPY
//

include { HAPPY_HAPPY      } from '../../../modules/nf-core/happy/happy/main'
include { HAPPY_PREPY      } from '../../../modules/nf-core/happy/prepy/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1    } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2    } from '../../../modules/local/bcftools/reheader'
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

    versions        = Channel.empty()
    tagged_variants = Channel.empty()

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

    if (params.preprocess.contains("prepy")){

        // apply prepy if required
        HAPPY_PREPY(
            input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
                [ meta, vcf, _regionsbed ]
            },
            fasta,
            fai
        )
        versions = versions.mix(HAPPY_PREPY.out.versions.first())

        test_ch = HAPPY_PREPY.out.preprocessed_vcf
    }

    // apply happy method for benchmarking
    HAPPY_HAPPY(
        test_ch.join(truth_ch, failOnDuplicate:true, failOnMismatch:true),
        fasta,
        fai,
        falsepositive_bed,
        stratification_tsv,
        stratification_bed
    )
    versions = versions.mix(HAPPY_HAPPY.out.versions.first())

    // tag meta and collect summary reports
    HAPPY_HAPPY.out.summary_csv
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "happy"], file) }
        .groupTuple()
        .set{ summary_reports }

    // Subsample TRUTH column from happy results
    BCFTOOLS_VIEW_TRUTH(
        HAPPY_HAPPY.out.vcf.join(HAPPY_HAPPY.out.tbi),
        [],
        [],
        []
        )
    versions = versions.mix(BCFTOOLS_VIEW_TRUTH.out.versions.first())

    // reheader benchmarking results properly and tag meta
    BCFTOOLS_REHEADER_1(
        BCFTOOLS_VIEW_TRUTH.out.vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_1.out.versions)

    // Subsample QUERY column from happy results
    BCFTOOLS_VIEW_QUERY(
        HAPPY_HAPPY.out.vcf.join(HAPPY_HAPPY.out.tbi),
        [],
        [],
        []
    )
    versions = versions.mix(BCFTOOLS_VIEW_QUERY.out.versions)

    // reheader benchmarking results properly and tag meta
    BCFTOOLS_REHEADER_2(
        BCFTOOLS_VIEW_QUERY.out.vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_2.out.versions)

    BCFTOOLS_FILTER_TRUTH_TP(
        BCFTOOLS_REHEADER_1.out.vcf.join(BCFTOOLS_REHEADER_1.out.index)
    )
    versions = versions.mix(BCFTOOLS_FILTER_TRUTH_TP.out.versions)

    BCFTOOLS_FILTER_TRUTH_TP.out.vcf
        .join(BCFTOOLS_FILTER_TRUTH_TP.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "happy"], file, index) }
        .set { vcf_tp_comp }

    BCFTOOLS_FILTER_TRUTH_FN(
        BCFTOOLS_REHEADER_1.out.vcf.join(BCFTOOLS_REHEADER_1.out.index)
    )
    versions = versions.mix(BCFTOOLS_FILTER_TRUTH_FN.out.versions)

    BCFTOOLS_FILTER_TRUTH_FN.out.vcf
        .join(BCFTOOLS_FILTER_TRUTH_FN.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "happy"], file, index) }
        .set { vcf_fn }

    BCFTOOLS_FILTER_QUERY_TP(
        BCFTOOLS_REHEADER_2.out.vcf.join(BCFTOOLS_REHEADER_2.out.index)
    )
    versions = versions.mix(BCFTOOLS_FILTER_QUERY_TP.out.versions)

    BCFTOOLS_FILTER_QUERY_TP.out.vcf
        .join(BCFTOOLS_FILTER_QUERY_TP.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "happy"], file, index) }
        .set { vcf_tp_base }

    BCFTOOLS_FILTER_QUERY_FP(
        BCFTOOLS_REHEADER_2.out.vcf.join(BCFTOOLS_REHEADER_2.out.index)
    )
    versions = versions.mix(BCFTOOLS_FILTER_QUERY_FP.out.versions)

    BCFTOOLS_FILTER_QUERY_FP.out.vcf
        .join(BCFTOOLS_FILTER_QUERY_FP.out.tbi)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "happy"], file, index) }
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
    versions        // channel: [versions.yml]

}
