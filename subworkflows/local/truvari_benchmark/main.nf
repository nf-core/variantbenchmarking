//
// TRUVARI_BENCHMARK: SUBWORKFLOW FOR TRUVARI BENCHMARKS
//

include { TRUVARI_BENCH                             } from '../../../modules/nf-core/truvari/bench'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1  } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2  } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_3  } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_4  } from '../../../modules/local/bcftools/reheader'

workflow TRUVARI_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    versions        = Channel.empty()
    tagged_variants = Channel.empty()

    TRUVARI_BENCH(
        input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
                [ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed ]
            },
        fasta,
        fai
    )
    versions = versions.mix(TRUVARI_BENCH.out.versions.first())

    TRUVARI_BENCH.out.summary
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "truvari"], file) }
        .groupTuple()
        .set { report }

    // reheader fn vcf files for tagged results
    BCFTOOLS_REHEADER_1(
        TRUVARI_BENCH.out.fn_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_1.out.versions)

    BCFTOOLS_REHEADER_1.out.vcf
        .join(BCFTOOLS_REHEADER_1.out.index)
        .map { _meta, file, _index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "truvari"], file) }
        .set { vcf_fn }

    // reheader fp vcf files for tagged results
    BCFTOOLS_REHEADER_2(
        TRUVARI_BENCH.out.fp_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_2.out.versions)

    // add tag and to meta
    BCFTOOLS_REHEADER_2.out.vcf
        .join(BCFTOOLS_REHEADER_2.out.index)
        .map { _meta, file, _index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "truvari"], file) }
        .set { vcf_fp }

    // reheader base tp vcf files for tagged results
    BCFTOOLS_REHEADER_3(
        TRUVARI_BENCH.out.tp_base_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_3.out.versions)

    // add tag and to meta
    BCFTOOLS_REHEADER_3.out.vcf
        .join(BCFTOOLS_REHEADER_3.out.index)
        .map { _meta, file, _index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "truvari"], file) }
        .set { vcf_tp_base }

    // reheader comp tp vcf files for tagged results
    BCFTOOLS_REHEADER_4(
        TRUVARI_BENCH.out.tp_comp_vcf.map{ meta, vcf ->
        [ meta, vcf, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_4.out.versions)

    // add tag and to meta
    BCFTOOLS_REHEADER_4.out.vcf
        .join(BCFTOOLS_REHEADER_4.out.index)
        .map { _meta, file, _index -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "truvari"], file) }
        .set { vcf_tp_comp }

    // collect tagged variant files
    tagged_variants = tagged_variants.mix(vcf_fn,
                                        vcf_fp,
                                        vcf_tp_base,
                                        vcf_tp_comp)

    emit:
    tagged_variants // channel: [val(meta), vcfs]
    report          // channel: [val(meta), reports]
    versions        // channel: [val(meta), versions.yml]
}
