//
// GERMLINE: SUBWORKFLOW FOR GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_BENCH       } from '../../modules/nf-core/truvari/bench'      addParams( options: params.options )

workflow GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index , truth_vcf, truth_index]
    bed       // channel: bed
    ref       // reference channel [ref.fa, ref.fa.fai],

    main:

    versions=Channel.empty()

    // SV benchmarking
    //
    // MODULE: TRUVARI_BENCH
    //
    //input_ch = input_ch.map{it -> tuple(it[0], it[1], it[2], it[3], it[4], [])}
    TRUVARI_BENCH(
        input_ch,
        bed,
        ref
    )
    versions = versions.mix(TRUVARI_BENCH.out.versions)


    emit:
    versions
}
