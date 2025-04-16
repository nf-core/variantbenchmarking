//
// SV_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SV GERMLINE VARIANTS
//

include { TRUVARI_BENCH            } from '../../../modules/nf-core/truvari/bench'
include { SVANALYZER_SVBENCHMARK   } from '../../../modules/nf-core/svanalyzer/svbenchmark'
include { WITTYER                  } from '../../../modules/nf-core/wittyer'
include { RTGTOOLS_BNDEVAL         } from '../../../modules/nf-core/rtgtools/bndeval'
include { RTGTOOLS_SVDECOMPOSE     } from '../../../modules/nf-core/rtgtools/svdecompose'
include { TABIX_BGZIP as TABIX_BGZIP_QUERY          } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_BGZIP as TABIX_BGZIP_TRUTH          } from '../../../modules/nf-core/tabix/bgzip'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1  } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2  } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_3  } from '../../../modules/local/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_4  } from '../../../modules/local/bcftools/reheader'

workflow SV_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    versions        = Channel.empty()
    summary_reports = Channel.empty()
    tagged_variants = Channel.empty()

    // SV benchmarking
    if (params.method.contains('truvari')){

        // use truvari benchmarking
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

        summary_reports = summary_reports.mix(report)

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
    }

    if (params.method.contains('svanalyzer') && params.variant_type != "copynumber"){
        // svbenchmark cannot be run with copynumber analysis
        // apply svanalyzer to benchmark SVs
        SVANALYZER_SVBENCHMARK(
            input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed, _targets_bed  ->
                    [ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _regionsbed ]
                },
            fasta,
            fai
        )
        versions = versions.mix(SVANALYZER_SVBENCHMARK.out.versions.first())

        // tag and collect summary file
        SVANALYZER_SVBENCHMARK.out.report
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "svbenchmark"], file) }
            .groupTuple()
            .set{ report }

        summary_reports = summary_reports.mix(report)

        // reheader fn vcf files for tagged results
        SVANALYZER_SVBENCHMARK.out.fns
            .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "svbenchmark"], file) }
            .set { vcf_fn }


        SVANALYZER_SVBENCHMARK.out.fps
            .map { _meta, file -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "svbenchmark"], file) }
            .set { vcf_fp }

        tagged_variants = tagged_variants.mix(
            vcf_fn,
            vcf_fp
        )

    }
    if (params.method.contains('wittyer')){

        // unzip vcf.gz files
        TABIX_BGZIP_QUERY(
            input_ch.map { meta, vcf, _tbi, _truth_vcf, _truth_tbi, _bed, _targets_bed  ->
                [ meta, vcf ]
            }
        )
        versions = versions.mix(TABIX_BGZIP_QUERY.out.versions.first())

        TABIX_BGZIP_TRUTH(
            input_ch.map { meta, _vcf, _tbi, truth_vcf, _truth_tbi, _bed, _targets_bed  ->
                [ meta, truth_vcf ]
            }
        )
        versions = versions.mix(TABIX_BGZIP_TRUTH.out.versions.first())

        input_ch.map { meta, _vcf, _tbi, _truth_vcf, _truth_tbi, bed, _targets_bed  ->
                [ meta, bed ]
            }
            .set { bed }

        //
        // MODULE: WITTYER
        //
        WITTYER(
            TABIX_BGZIP_QUERY.out.output
                .join(TABIX_BGZIP_TRUTH.out.output, failOnDuplicate:true, failOnMismatch:true)
                .join(bed, failOnDuplicate:true, failOnMismatch:true)
        )
        versions = versions.mix(WITTYER.out.versions.first())

        WITTYER.out.report
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "wittyer"], file) }
            .groupTuple()
            .set{ report }
        summary_reports = summary_reports.mix(report)

    }
    if (params.method.contains('bndeval')) {

        RTGTOOLS_SVDECOMPOSE(
            input_ch.map { meta, _vcf, _tbi, truth_vcf, _truth_tbi, _bed, _targets_bed  ->
                [ meta, truth_vcf, _truth_tbi  ]
            }
        )

        input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _bed, _targets_bed  ->
                    [ meta, vcf, _tbi, _bed ]
                }.join(RTGTOOLS_SVDECOMPOSE.out.bnd)
                .map{meta, vcf, tbi, bed, bnd -> [meta, vcf, tbi, bnd, [], bed]}
                .set{bndeval_input}

        RTGTOOLS_BNDEVAL(
            bndeval_input
        )
    }

    emit:
    tagged_variants // channel: [val(meta), vcfs]
    summary_reports // channel: [val(meta), reports]
    versions        // channel: [val(meta), versions.yml]
}
