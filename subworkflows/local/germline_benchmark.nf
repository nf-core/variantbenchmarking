//
// GERMLINE: SUBWORKFLOW FOR GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_BENCH       } from '../../modules/nf-core/truvari/bench'      addParams( options: params.options )

workflow GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf,test_index, truth_vcf, truth_index, bed]
    ref       // reference channel [ref.fa, ref.fa.fai],

    main:

    versions=Channel.empty()

    // SV benchmarking
    //
    // MODULE: TRUVARI_BENCH
    //
    TRUVARI_BENCH(
        input_ch,
        ref
    )
    versions = versions.mix(TRUVARI_BENCH.out.versions)


    emit:
    versions
}
