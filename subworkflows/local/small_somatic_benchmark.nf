//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

params.options = [:]

include { HAPPY_SOMPY      } from '../../modules/nf-core/happy/sompy/main'               addParams( options: params.options )
include { HAPPY_FTX        } from '../../modules/local/happy_ftx'                  addParams( options: params.options )

workflow SMALL_SOMATIC_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf,test_index, truth_vcf, truth_index, bed]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions=Channel.empty()
    summary_reports=Channel.empty()

    if (params.method.contains('sompy')){
        // Sompy analysis requires feature table: https://github.com/Illumina/hap.py/blob/master/example/sompy/reference_outputs/sompy.strelka_grch38_admix_pass_indels.features.csv
        // I dont know how to prepare it for now.

        input_ch.map { it -> tuple(it[0], it[3], it[5], [], "generic") }
                .set{ch_ftx}

        HAPPY_FTX(
            ch_ftx,
            fasta,
            fai,
            [[],[]]
        )

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

    emit:
    summary_reports
    versions
}
