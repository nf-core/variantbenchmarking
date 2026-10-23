//
// SVANALYZER_BENCHMARK: SUBWORKFLOW FOR SVANALYZER BENCHMARK
//

include { SVANALYZER_SVBENCHMARK  } from '../../../modules/nf-core/svanalyzer/svbenchmark'
include { SUBTRACT_VCF as SUBTRACT_VCF_TRUTH  } from '../../../modules/local/custom/subtract_vcf/'
include { SUBTRACT_VCF as SUBTRACT_VCF_QUERY  } from '../../../modules/local/custom/subtract_vcf/'

workflow SVANALYZER_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    tagged_variants = channel.empty()

    // apply svanalyzer to benchmark SVs
    SVANALYZER_SVBENCHMARK(
        input_ch.map{ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed, _targets_bed  ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi, regionsbed ]
            },
        fasta,
        fai
    )

    // tag and collect summary file
    SVANALYZER_SVBENCHMARK.out.report
        .map { _meta, report -> tuple([vartype: params.variant_type] + [benchmark_tool: "svbenchmark"], report) }
        .groupTuple()
        .set{ report }

    // reheader fn vcf files for tagged results
    SVANALYZER_SVBENCHMARK.out.fns
        .map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "svbenchmark"], vcf) }
        .set { vcf_fn }

    SVANALYZER_SVBENCHMARK.out.fps
        .map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "svbenchmark"], vcf) }
        .set { vcf_fp }

    // subtract FPs from Query to find TPs in Query
    SUBTRACT_VCF_QUERY(
        input_ch.map{ meta, vcf, tbi, _truth_vcf, _truth_tbi, regionsbed, targets_bed  ->
            [ meta, vcf, tbi, regionsbed, targets_bed ]}.join(SVANALYZER_SVBENCHMARK.out.fps)
        )

    // subtract Fns from Truth to find TPs in TRUTH
    SUBTRACT_VCF_TRUTH(
        input_ch.map{ meta, _vcf, _tbi, truth_vcf, truth_tbi, regionsbed, targets_bed  ->
            [ meta, truth_vcf, truth_tbi, regionsbed, targets_bed ]}.join(SVANALYZER_SVBENCHMARK.out.fns)
        )

    // reheader tp_comp vcf files for tagged results
    SUBTRACT_VCF_QUERY.out.vcf
        .map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "svbenchmark"], vcf) }
        .set { vcf_tp_comp }

    // reheader tp_base vcf files for tagged results
    SUBTRACT_VCF_TRUTH.out.vcf
        .map { _meta, vcf -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "svbenchmark"], vcf) }
        .set { vcf_tp_base }

    tagged_variants = tagged_variants.mix(
        vcf_fn,
        vcf_fp,
        vcf_tp_comp,
        vcf_tp_base
    )

    emit:
    tagged_variants // channel: [val(meta), vcfs]
    report          // channel: [val(meta), reports]
}
