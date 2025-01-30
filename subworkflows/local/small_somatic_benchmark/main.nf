//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

include { HAPPY_SOMPY           } from '../../../modules/nf-core/happy/sompy/main'
include { BAMSURGEON_EVALUATOR  } from '../../../modules/local/bamsurgeon/evaluator'

workflow SMALL_SOMATIC_BENCHMARK {
    take:
    input_ch    // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, bed]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions        = Channel.empty()
    summary_reports = Channel.empty()

    if (params.method.contains('sompy')){
        // apply sompy for small somatic variant benchmarking
        HAPPY_SOMPY(
            input_ch.map { meta, vcf, _tbi, truth_vcf, _truth_tbi, bed ->
                [ meta, vcf, truth_vcf, bed, [] ]
            },
            fasta,
            fai,
            [[],[]],
            [[],[]],
            [[],[]]
        )
        versions = versions.mix(HAPPY_SOMPY.out.versions.first())

        HAPPY_SOMPY.out.stats
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "sompy"], file) }
            .groupTuple()
            .set{ report }
        summary_reports = summary_reports.mix(report)
    }

    // not working for now
    if (params.method.contains('bamsurgeon')){

        BAMSURGEON_EVALUATOR(
            input_ch.map { meta, vcf, tbi, truth_vcf, truth_tbi, _bed ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi ]
            },
            fasta,
            fai
        )
        versions = versions.mix(BAMSURGEON_EVALUATOR.out.versions.first())

        BAMSURGEON_EVALUATOR.out.stats
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "bamsurgeon"], file) }
            .groupTuple()
            .set{ report }
        summary_reports = summary_reports.mix(report)
    }

    emit:
    summary_reports // channel: [val(meta), reports]
    versions        // channel: [versions.yml]
}
