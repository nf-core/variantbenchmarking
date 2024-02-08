//
// GERMLINE: SUBWORKFLOW FOR GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_PHAB           } from '../../modules/local/truvari_phab'                  addParams( options: params.options )
include { TRUVARI_BENCH          } from '../../modules/nf-core/truvari/bench'          addParams( options: params.options )
include { SVANALYZER_SVBENCHMARK } from '../../modules/nf-core/svanalyzer/svbenchmark' addParams( options: params.options )
include { WITTYER                } from '../../modules/local/wittyer'                  addParams( options: params.options )
include { VCFDIST                } from '../../modules/local/vcfdist'                  addParams( options: params.options )
include { BAMSURGEON_EVALUATOR   } from '../../modules/local/bamsurgeon_evaluator'     addParams( options: params.options )

workflow GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),val(meta2), test_vcf, test_index , truth_vcf, truth_index]
    ref       // reference channel [ref.fa, ref.fa.fai]
    truth_vcf // channel: [val(meta),val(meta2),truth_vcf, truth_index]

    main:

    versions=Channel.empty()

    // SV benchmarking

    if (params.method.contains('truvari')){

        if(params.harmonize){
            //
            // TRUVARI: TRUVARI_PHAB
            //
            TRUVARI_PHAB(
                input_ch,
                ref
            )
        }
        //
        // MODULE: TRUVARI_BENCH
        //
        TRUVARI_BENCH(
            input_ch,
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
            ref
            )
        versions = versions.mix(SVANALYZER_SVBENCHMARK.out.versions)

    } 

    if (params.method.contains('wittyer')){
        //
        // MODULE: WITTYER
        //
        WITTYER(
            input_ch,
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
            ref
        )
        versions = versions.mix(VCFDIST.out.versions)
    }

    if (params.method.contains('bamsurgeon')){
        //
        // MODULE: BAMSURGEON_EVALUATOR
        //
        //https://github.com/adamewing/bamsurgeon/blob/master/scripts/evaluator.py
        BAMSURGEON_EVALUATOR(
            input_ch.map{it -> tuple(it[0],it[1], it[2], it[3], it[4], it[5])},
            ref,
            "SV"
        )
        versions = versions.mix(BAMSURGEON_EVALUATOR.out.versions)
    }



    emit:
    versions
}
