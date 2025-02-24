//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

include { HAPPY_SOMPY           } from '../../../modules/nf-core/happy/sompy/main'

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
    emit:
    summary_reports // channel: [val(meta), reports]
    versions        // channel: [versions.yml]
}
