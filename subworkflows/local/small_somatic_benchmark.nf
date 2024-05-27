//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

params.options = [:]

include { HAPPY_SOMPY      } from '../../modules/nf-core/happy/sompy/main'               addParams( options: params.options )
include { BAMSURGEON_EVALUATOR  } from '../../modules/local/bamsurgeon_evaluator.nf'     addParams( options: params.options )

workflow SMALL_SOMATIC_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, bed]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions=Channel.empty()
    summary_reports=Channel.empty()

    if (params.method.contains('sompy')){

        HAPPY_SOMPY(
            input_ch.map { it -> tuple(it[0], it[3], it[1], it[5], []) },
            fasta,
            fai,
            [[],[]],
            [[],[]],
            [[],[]]
        )
        versions = versions.mix(HAPPY_SOMPY.out.versions)

        HAPPY_SOMPY.out.stats
            .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "sompy"], file) }
            .groupTuple()
            .set{ report}
        summary_reports = summary_reports.mix(report)
    }

    if (params.method.contains('bamsurgeon')){

        BAMSURGEON_EVALUATOR(
            input_ch.map { it -> tuple(it[0], it[1], it[2], it[3], it[4]) },
            fasta,
            fai
        )
        versions = versions.mix(BAMSURGEON_EVALUATOR.out.versions)

        BAMSURGEON_EVALUATOR.out.stats
            .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "bamsurgeon"], file) }
            .groupTuple()
            .set{ report}
        summary_reports = summary_reports.mix(report)
    }

    emit:
    summary_reports
    versions
}
