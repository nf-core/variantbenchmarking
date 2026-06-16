//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

include { SOMPY_BENCHMARK     } from '../../../subworkflows/local/sompy_benchmark'
include { RTGTOOLS_VCFEVAL    } from '../../../modules/nf-core/rtgtools/vcfeval'


workflow SMALL_SOMATIC_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index,  regionsbed, targets_bed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    sdf                // reference channel [val(meta), sdf]
    falsepositive_bed  // reference channel [val(meta), bed]
    ambiguous_beds     // reference channel [val(meta), bed]

    main:

    summary_reports     = channel.empty()
    tagged_variants     = channel.empty()
    tagged_variants_csv = channel.empty()

    if (params.method.contains('sompy')){
        SOMPY_BENCHMARK(
            input_ch,
            fasta,
            fai,
            falsepositive_bed,
            ambiguous_beds
        )
        summary_reports     = summary_reports.mix(SOMPY_BENCHMARK.out.summary_reports)
        tagged_variants_csv = tagged_variants_csv.mix(SOMPY_BENCHMARK.out.tagged_variants_csv)
                                    .map { meta, file -> tuple(meta + [vartype: params.variant_type] + [id: "sompy"], file) }
    }

    if (params.method.contains('rtgtools')){

        RTGTOOLS_VCFEVAL(
            input_ch,
            sdf
        )
        summary_reports = summary_reports.mix(RTGTOOLS_VCFEVAL.out.summary
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], file) }
        .groupTuple())
        tagged_variants = tagged_variants.mix(RTGTOOLS_VCFEVAL.out.fn_vcf.join(RTGTOOLS_VCFEVAL.out.fn_tbi),
                                            RTGTOOLS_VCFEVAL.out.fp_vcf.join(RTGTOOLS_VCFEVAL.out.fp_tbi),
                                            RTGTOOLS_VCFEVAL.out.baseline_vcf.join(RTGTOOLS_VCFEVAL.out.baseline_tbi),
                                            RTGTOOLS_VCFEVAL.out.tp_vcf.join(RTGTOOLS_VCFEVAL.out.tp_tbi))
                                        .map { _meta, file, index ->
                                            def mapping = [
                                                'fn': 'FN',
                                                'fp': 'FP',
                                                'tp-baseline': 'TP_base',
                                                'tp': 'TP_comp'
                                            ]
                                            def tag = file.getName().tokenize('.').find { token -> token in ['fn', 'fp', 'tp-baseline', 'tp'] }
                                            def transformedTag = mapping[tag] ?: tag
                                            tuple([vartype: params.variant_type, id: "rtgtools", tag: transformedTag], file, index)
                                        }
    }

    emit:
    summary_reports     // channel: [val(meta), reports]
    tagged_variants     // channel: [val(meta), vcfs]
    tagged_variants_csv // channel: [val(meta), csvs]
}
