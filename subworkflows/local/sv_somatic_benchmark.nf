//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

include { BAMSURGEON_EVALUATOR  } from '../../modules/local/bamsurgeon_evaluator.nf'

workflow SV_SOMATIC_BENCHMARK {
    take:
    input_ch    // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, bed]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions =          Channel.empty()
    summary_reports =   Channel.empty()


    if (params.method.contains('bamsurgeon')){

        BAMSURGEON_EVALUATOR(
            input_ch.map { meta, vcf, tbi, truth_vcf, truth_tbi, bed ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi ]
            },
            fasta,
            fai
        )
        versions = versions.mix(BAMSURGEON_EVALUATOR.out.versions.first())

        BAMSURGEON_EVALUATOR.out.stats
            .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "bamsurgeon"], file) }
            .groupTuple()
            .set{ report }
        summary_reports = summary_reports.mix(report)
    }

    emit:
    summary_reports
    versions
}
