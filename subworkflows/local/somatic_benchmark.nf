//
// SOMATIC: SUBWORKFLOW FOR SOMATIC VARIANTS
//

params.options = [:]

include { TRUVARI_BENCH                       } from '../../modules/nf-core/truvari/bench'          addParams( options: params.options )
include { HAPPY_SOMPY                         } from '../../modules/nf-core/happy/sompy'            addParams( options: params.options )
include { SVANALYZER_SVBENCHMARK              } from '../../modules/nf-core/svanalyzer/svbenchmark' addParams( options: params.options )

workflow SOMATIC_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf,test_index, truth_vcf, truth_index]
    ref       // reference channel [ref.fa, ref.fa.fai]
    sv_bed
    snv_bed

    main:

    versions=Channel.empty()

    input_ch.branch{
        sv: it[0].vartype == "sv"
        indel: it[0].vartype == "indel"
        snv: it[0].vartype == "snv"
        }.set{bench}

    // SV Benchmarking
    //
    // MODULE: TRUVARI_BENCH
    //
    TRUVARI_BENCH(
        bench.sv,
        sv_bed,
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

    HAPPY_SOMPY(
        bench.snv,
        snv_bed,
        [],
        ref,
        [[],[]],
        [[],[]],
        [[],[]]
    )
    versions = versions.mix(HAPPY_SOMPY.out.versions)

    emit:
    versions
}
