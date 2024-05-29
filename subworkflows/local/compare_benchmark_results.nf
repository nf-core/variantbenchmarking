
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge tp/fp/fn results from different tools.
//

params.options = [:]

include { SURVIVOR_MERGE  } from '../../modules/nf-core/survivor/merge'  addParams( options: params.options )
include { TABIX_BGZIP     } from '../../modules/nf-core/tabix/bgzip'     addParams( options: params.options )

workflow COMPARE_BENCHMARK_RESULTS {
    take:
    input_ch    // channel: [val(meta), vcf.gz]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:
    versions   = Channel.empty()

    //
    // MODULE: TABIX_BGZIP
    //
    // unzip

    TABIX_BGZIP(
        input_ch
    )


    emit:
    versions
}
