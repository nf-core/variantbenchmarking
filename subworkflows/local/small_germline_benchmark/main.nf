//
// SMALL_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SMALL GERMLINE VARIANTS
//

include { RTGTOOLS_FORMAT  } from '../../../modules/nf-core/rtgtools/format/main'
include { RTGTOOLS_VCFEVAL } from '../../../modules/nf-core/rtgtools/vcfeval/main'
include { HAPPY_HAPPY      } from '../../../modules/nf-core/happy/happy/main'
include { HAPPY_PREPY      } from '../../../modules/nf-core/happy/prepy/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_1  } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_2  } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_3  } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_4  } from '../../../modules/nf-core/bcftools/reheader'

workflow SMALL_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),test_vcf,test_index,truth_vcf,truth_index, bed]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]
    sdf       // reference channel [val(meta), sdf]

    main:

    versions        = Channel.empty()
    summary_reports = Channel.empty()
    tagged_variants = Channel.empty()

    if (params.method.contains('rtgtools')){

        if (!params.sdf){

            // Use rtgtools format to generate sdf file if necessary
            RTGTOOLS_FORMAT(
                fasta.map { meta, file -> [ meta, file, [], [] ] }
            )
            versions = versions.mix(RTGTOOLS_FORMAT.out.versions)
            sdf = RTGTOOLS_FORMAT.out.sdf

        }

        // apply rtgtools eval method
        RTGTOOLS_VCFEVAL(
            input_ch.map { meta, vcf, tbi, truth_vcf, truth_tbi, bed ->
                [ meta, vcf, tbi, truth_vcf, truth_tbi, bed, [] ]
            },
            sdf
        )
        versions = versions.mix(RTGTOOLS_VCFEVAL.out.versions.first())

        // collect summary reports
        RTGTOOLS_VCFEVAL.out.summary
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "rtgtools"], file) }
            .groupTuple()
            .set{ report }

        summary_reports = summary_reports.mix(report)

        // reheader benchmarking results properly and tag meta
        BCFTOOLS_REHEADER_1(
            RTGTOOLS_VCFEVAL.out.fn_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_1.out.versions.first())

        BCFTOOLS_REHEADER_1.out.vcf
            .join(BCFTOOLS_REHEADER_1.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "rtgtools"], file, index) }
            .set { vcf_fn }

        BCFTOOLS_REHEADER_2(
            RTGTOOLS_VCFEVAL.out.fp_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_2.out.versions)

        BCFTOOLS_REHEADER_2.out.vcf
            .join(BCFTOOLS_REHEADER_2.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "rtgtools"], file, index) }
            .set { vcf_fp }

        BCFTOOLS_REHEADER_3(
            RTGTOOLS_VCFEVAL.out.baseline_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_3.out.versions)

        BCFTOOLS_REHEADER_3.out.vcf
            .join(BCFTOOLS_REHEADER_3.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "rtgtools"], file, index) }
            .set { vcf_tp_base }

        BCFTOOLS_REHEADER_4(
            RTGTOOLS_VCFEVAL.out.tp_vcf.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
            },
            fai
        )
        versions = versions.mix(BCFTOOLS_REHEADER_4.out.versions)

        BCFTOOLS_REHEADER_4.out.vcf
            .join(BCFTOOLS_REHEADER_4.out.index)
            .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "rtgtools"], file, index) }
            .set { vcf_tp_comp }

        tagged_variants = tagged_variants.mix(
            vcf_fn,
            vcf_fp,
            vcf_tp_base,
            vcf_tp_comp
        )
    }

    if (params.method.contains('happy')){

        input_ch
            .map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _bed ->
                [ meta, vcf ]
            }
            .set { test_ch }

        input_ch
            .map{ meta, _vcf, _tbi, truth_vcf, _truth_tbi, bed ->
                [ meta, truth_vcf, bed, [] ]
            }
            .set { truth_ch }

        if (params.preprocess.contains("prepy")){

            // apply prepy if required
            HAPPY_PREPY(
                input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, bed ->
                    [ meta, vcf, bed ]
                },
                fasta,
                fai
            )
            versions = versions.mix(HAPPY_PREPY.out.versions.first())
            // TODO: Check norm settings https://github.com/Illumina/hap.py/blob/master/doc/normalisation.md

            test_ch = HAPPY_PREPY.out.preprocessed_vcf
        }

        // apply happy method for benchmarking
        HAPPY_HAPPY(
            test_ch.join(truth_ch, failOnDuplicate:true, failOnMismatch:true),
            fasta,
            fai,
            [[],[]],
            [[],[]],
            [[],[]]
        )
        versions = versions.mix(HAPPY_HAPPY.out.versions.first())

        // tag meta and collect summary reports
        HAPPY_HAPPY.out.summary_csv
            .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "happy"], file) }
            .groupTuple()
            .set{ report }
        summary_reports = summary_reports.mix(report)
    }
    emit:
    summary_reports // channel: [val(meta), reports]
    tagged_variants // channel: [val(meta), vcfs]
    versions        // channel: [versions.yml]

}
