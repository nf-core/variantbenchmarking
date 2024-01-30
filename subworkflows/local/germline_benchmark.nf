//
// GERMLINE: SUBWORKFLOW FOR GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_BENCH          } from '../../modules/nf-core/truvari/bench'          addParams( options: params.options )
include { SVANALYZER_SVBENCHMARK } from '../../modules/nf-core/svanalyzer/svbenchmark' addParams( options: params.options )
include { WITTYER                } from '../../modules/local/wittyer'                  addParams( options: params.options )
include { VCFDIST                } from '../../modules/local/vcfdist'                  addParams( options: params.options )

workflow GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index , truth_vcf, truth_index]
    bed       // channel: bed
    ref       // reference channel [ref.fa, ref.fa.fai],

    main:

    versions=Channel.empty()

    // SV benchmarking

    if (params.method.contains('truvari')){
        //
        // MODULE: TRUVARI_BENCH
        //
        TRUVARI_BENCH(
            input_ch,
            bed,
            ref
        )
        versions = versions.mix(TRUVARI_BENCH.out.versions)
    }

    if (params.method.contains('svanalyzer')){
        //
        // MODULE: SVANALYZER_SVBENCHMARK
        //
        // slower than truvari
        SVANALYZER_SVBENCHMARK(
            input_ch,
            ref,
            bed        
        )
        versions = versions.mix(SVANALYZER_SVBENCHMARK.out.versions)
    } 

    if (params.method.contains('wittyer')){
        //
        // MODULE: WITTYER
        //
        WITTYER(
            input_ch,
            bed,
            []
        )
        versions = versions.mix(WITTYER.out.versions)
    }


    if (params.method.contains('vcfdist')){
        //
        // MODULE: VCFDIST
        //
        VCFDIST(
            input_ch,
            ref,
            bed
        )
        versions = versions.mix(VCFDIST.out.versions)
    }


    emit:
    versions
}
