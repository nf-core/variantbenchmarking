//
// GERMLINE: SUBWORKFLOW FOR GERMLINE VARIANTS
//

params.options = [:]

include { PREPARE_VCFS as PREPARE_VCFS_TRUTH  } from '../../modules/local/prepare_vcfs.nf'     addParams( options: params.options )
include { PREPARE_VCFS as PREPARE_VCFS_TEST   } from '../../modules/local/prepare_vcfs.nf'     addParams( options: params.options )
include { TRUVARI_BENCH                       } from '../../modules/nf-core/truvari/bench'     addParams( options: params.options )
include { HAPPY_HAPPY                         } from '../../modules/nf-core/happy/happy'        addParams( options: params.options )

workflow GERMLINE {
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

    //
    // MODULE: TRUVARI_BENCH
    //
    TRUVARI_BENCH(
        bench.sv,
        sv_bed,
        ref
    )
    versions = versions.mix(TRUVARI_BENCH.out.versions)

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
