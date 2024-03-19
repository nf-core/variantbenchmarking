//
// SMALL_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SMALL GERMLINE VARIANTS
//

params.options = [:]

include { RTGTOOLS_FORMAT  } from '../../modules/nf-core/rtgtools/format/main'           addParams( options: params.options )
include { RTGTOOLS_VCFEVAL } from '../../modules/nf-core/rtgtools/vcfeval/main'          addParams( options: params.options )
include { HAPPY_HAPPY      } from '../../modules/nf-core/happy/happy/main'               addParams( options: params.options )

workflow SMALL_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),test_vcf,test_index,truth_vcf,truth_index, bed]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions=Channel.empty()

    if (params.method.contains('rtgtools')){
        //
        // MODULE: RTGTOOLS_FORMAT
        //
        RTGTOOLS_FORMAT(
            fasta.map { it -> tuple([id: it[1].getSimpleName(), "pair": "single_end"], it[1], [], []) }
        )
        versions = versions.mix(RTGTOOLS_FORMAT.out.versions)

        //
        // MODULE: RTGTOOLS_VCFEVAL
        //
        RTGTOOLS_VCFEVAL(
            input_ch.map { it -> tuple(it[0],it[3], it[4], it[1], it[2], it[5], []) },
            RTGTOOLS_FORMAT.out.sdf
        )
        versions = versions.mix(RTGTOOLS_VCFEVAL.out.versions)
    }

    if (params.method.contains('happy')){

        HAPPY_HAPPY(
            input_ch.map { it -> tuple(it[0],it[3], it[1], it[5], []) },
            fasta,
            fai,
            [[],[]],
            [[],[]],
            [[],[]]
        )
        versions = versions.mix(HAPPY_HAPPY.out.versions)
    }

    emit:
    versions
}
