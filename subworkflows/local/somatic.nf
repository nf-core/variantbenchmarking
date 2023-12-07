//
// SOMATIC: SUBWORKFLOW FOR SOMATIC VARIANTS
//

params.options = [:]

include { PREPARE_VCFS as PREPARE_VCFS_TRUTH  } from '../../modules/local/prepare_vcfs.nf'     addParams( options: params.options )
include { PREPARE_VCFS as PREPARE_VCFS_TEST   } from '../../modules/local/prepare_vcfs.nf'     addParams( options: params.options )
include { TRUVARI_BENCH                       } from '../../modules/nf-core/truvari/bench'     addParams( options: params.options )
include { HAPPY_SOMPY                         } from '../../modules/nf-core/happy/sompy'     addParams( options: params.options )

workflow SOMATIC {
    take:
    input_ch  // channel: [val(meta), test_vcf, truth_vcf]
    ref       // reference channel [ref.fa, ref.fa.fai]
    sv_bed
    snv_bed

    main:

    versions=Channel.empty()

    //
    // PREPARE_VCFS
    //
    input_ch.map{ it -> [it[0], it[2]]}
            .set{truth_vcf}

    input_ch.map{ it -> [it[0], it[1]]}
            .set{test_vcf}

    PREPARE_VCFS_TRUTH(
        truth_vcf
    )

    PREPARE_VCFS_TEST(
        test_vcf
    )
    versions = versions.mix(PREPARE_VCFS_TEST.out.versions)

    bench_ch = PREPARE_VCFS_TEST.out.gz_tbi.join(PREPARE_VCFS_TRUTH.out.gz_tbi)

    bench_ch.branch{
        sv: it[0].vartype == "sv"
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
