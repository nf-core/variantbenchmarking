//
// SOMATIC: SUBWORKFLOW FOR SOMATIC VARIANTS
//

params.options = [:]

include { TRUVARI_BENCH                       } from '../../modules/nf-core/truvari/bench'          addParams( options: params.options )
include { SVANALYZER_SVBENCHMARK              } from '../../modules/nf-core/svanalyzer/svbenchmark' addParams( options: params.options )

workflow SOMATIC_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf,test_index, truth_vcf, truth_index, bed]
    ref       // reference channel [ref.fa, ref.fa.fai]

    main:

    versions=Channel.empty()

    // SV Benchmarking
    //
    // MODULE: TRUVARI_BENCH
    //
    TRUVARI_BENCH(
        input_ch,
        ref
    )
    versions = versions.mix(TRUVARI_BENCH.out.versions)

    // SV Benchmarking
    //
    // MODULE: SVANALYZER_SVBENCHMARK
    //
    // note: slow
    //SVANALYZER_SVBENCHMARK(
    //    bench.sv,
    //    ref,
    //    sv_bed
    //)
    //versions = versions.mix(SVANALYZER_SVBENCHMARK.out.versions)

    // Small Variant Benchmarking

    // SOMPY https://sites.google.com/view/seqc2/home/benchmarking-examples?authuser=0
    // is used for somatic variant benchmarking!


    emit:
    versions
}
