//
// GERMLINE: SUBWORKFLOW FOR GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_BENCH       } from '../../modules/nf-core/truvari/bench'      addParams( options: params.options )
include { HAPPY_HAPPY         } from '../../modules/nf-core/happy/happy'        addParams( options: params.options )
include { RTGTOOLS_VCFEVAL    } from '../../modules/nf-core/rtgtools/vcfeval'   addParams( options: params.options )
include { RTGTOOLS_FORMAT     } from '../../modules/nf-core/rtgtools/format'    addParams( options: params.options )
include { RTGTOOLS_ROCPLOT    } from '../../modules/nf-core/rtgtools/rocplot'   addParams( options: params.options )

workflow GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf,test_index, truth_vcf, truth_index]
    ref       // reference channel [ref.fa, ref.fa.fai],
    sdf_file
    sv_bed
    snv_bed

    main:

    versions=Channel.empty()

    input_ch.branch{
        sv: it[0].vartype == "sv"
        indel: it[0].vartype == "indel"
        snv: it[0].vartype == "snv"
        }.set{bench}

    // SV benchmarking
    //
    // MODULE: TRUVARI_BENCH
    //
    TRUVARI_BENCH(
        bench.sv,
        sv_bed,
        ref
    )
    versions = versions.mix(TRUVARI_BENCH.out.versions)

    // Small variant benchmarking


    //
    // MODULE: RTGTOOLS_VCFEVAL
    //
    if (!params.sdf_file){
        ref.map { it -> tuple([id: it[0].baseName, single_end: "true"], it[0], []) }
            .set{fasta}

        RTGTOOLS_FORMAT(
            fasta,
            []
        )
        versions = versions.mix(RTGTOOLS_FORMAT.out.versions)
        sdf_file = RTGTOOLS_FORMAT.out.sdf
    }

    RTGTOOLS_VCFEVAL(
        bench.snv,
        snv_bed,
        [],
        sdf_file
    )
    versions = versions.mix(RTGTOOLS_VCFEVAL.out.versions)

    RTGTOOLS_ROCPLOT(
        RTGTOOLS_VCFEVAL.out.weighted_roc
    )
    versions = versions.mix(RTGTOOLS_ROCPLOT.out.versions)

    //
    // MODULE: HAPPY
    //
    HAPPY_HAPPY(
        bench.snv,
        snv_bed,
        [],
        ref,
        [[],[]],
        [[],[]],
        [[],[]]
    )
    versions = versions.mix(HAPPY_HAPPY.out.versions)

    emit:
    versions
}
