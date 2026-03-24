//
// SOMATIC: SUBWORKFLOW FOR SOMPY BENCHMARKING
//

include { HAPPY_SOMPY          } from '../../../modules/nf-core/happy/sompy'
include { SOMPY_FEATURES_SPLIT } from '../../../modules/local/sompy_features/split'

workflow SOMPY_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index,  regionsbed, targets_bed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    falsepositive_bed  // reference channel [val(meta), bed]
    ambiguous_beds     // reference channel [val(meta), bed]

    main:

    tagged_variants_csv = channel.empty()

    // apply sompy for small somatic variant benchmarking
    HAPPY_SOMPY(
        input_ch.map{meta, test, test_index, truth, truth_index, regions, target -> [meta, test, truth, regions, target]},
        fasta,
        fai,
        falsepositive_bed,
        ambiguous_beds,
        [[],[]]
    )

    HAPPY_SOMPY.out.stats
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "sompy"], file) }
        .groupTuple()
        .set{ summary_reports }

    SOMPY_FEATURES_SPLIT(
        HAPPY_SOMPY.out.features
    )

    SOMPY_FEATURES_SPLIT.out.TP
        .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "sompy"], file) }
        .set{tp_vars}

    SOMPY_FEATURES_SPLIT.out.FP
        .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "sompy"], file) }
        .set{fp_vars}

    SOMPY_FEATURES_SPLIT.out.FN
        .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "sompy"], file) }
        .set{fn_vars}

    tagged_variants_csv = tagged_variants_csv
                            .mix(tp_vars)
                            .mix(fp_vars)
                            .mix(fn_vars)

    emit:
    summary_reports     // channel: [val(meta), reports]
    tagged_variants_csv // channel: [val(meta), csvs]
}
