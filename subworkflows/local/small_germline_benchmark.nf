//
// SMALL_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SMALL GERMLINE VARIANTS
//

params.options = [:]

include { RTGTOOLS_FORMAT  } from '../../modules/nf-core/rtgtools/format/main'           addParams( options: params.options )
include { RTGTOOLS_VCFEVAL } from '../../modules/nf-core/rtgtools/vcfeval/main'          addParams( options: params.options )
include { HAPPY_HAPPY      } from '../../modules/nf-core/happy/happy/main'               addParams( options: params.options )
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_SUBSET } from '../../modules/nf-core/bcftools/view/main' addParams( options: params.options )

workflow SMALL_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),test_vcf,test_index,truth_vcf,truth_index, bed]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]
    sdf       // reference channel [val(meta), sdf]

    main:

    versions        =Channel.empty()
    summary_reports =Channel.empty()
    tagged_variants =Channel.empty()

    if (params.method.contains('rtgtools')){

        if (!params.sdf){
            //
            // MODULE: RTGTOOLS_FORMAT
            //
            RTGTOOLS_FORMAT(
                fasta.map { it -> tuple([id: it[1].getSimpleName()], it[1], [], []) }
            )
            versions = versions.mix(RTGTOOLS_FORMAT.out.versions)
            sdf = RTGTOOLS_FORMAT.out.sdf
        }
        test_ch = input_ch.map { it -> tuple(it[0], it[1], [2]) }
        if (params.preprocess.contains("prepy")){
            //
            // MODULE: BCFTOOLS_VIEW
            //
            BCFTOOLS_VIEW_SUBSET(
                input_ch.map { it -> tuple(it[0], it[1], it[2]) },
                [],[],[]
            )
            versions = versions.mix(BCFTOOLS_VIEW_SUBSET.out.versions)

            BCFTOOLS_VIEW_SUBSET.out.vcf.map{it -> tuple(it[0], it[1], [])}.set{test_ch}
        }

        truth_ch = input_ch.map { it -> tuple(it[0], it[3], it[4], it[5], []) }

        //
        // MODULE: RTGTOOLS_VCFEVAL
        //
        RTGTOOLS_VCFEVAL(
            test_ch.join(truth_ch),
            sdf
        )
        versions = versions.mix(RTGTOOLS_VCFEVAL.out.versions)

        RTGTOOLS_VCFEVAL.out.summary
            .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "rtgtools"], file) }
            .groupTuple()
            .set{ report}

        summary_reports = summary_reports.mix(report)

        RTGTOOLS_VCFEVAL.out.fn_vcf
            .join(RTGTOOLS_VCFEVAL.out.fn_tbi)
            .map { meta, file, index -> tuple([vartype: meta.vartype] + [tag: "FN"] + [id: "rtgtools"], file, index) }
            .set { vcf_fn }

        RTGTOOLS_VCFEVAL.out.fp_vcf
            .join(RTGTOOLS_VCFEVAL.out.fp_tbi)
            .map { meta, file, index -> tuple([vartype: meta.vartype] + [tag: "FP"] + [id: "rtgtools"], file, index) }
            .set { vcf_fp }

        RTGTOOLS_VCFEVAL.out.baseline_vcf
            .join(RTGTOOLS_VCFEVAL.out.baseline_tbi)
            .map { meta, file, index -> tuple([vartype: meta.vartype] + [tag: "TP_base"] + [id: "rtgtools"], file, index) }
            .set { vcf_tp_base }

        RTGTOOLS_VCFEVAL.out.tp_vcf
            .join(RTGTOOLS_VCFEVAL.out.tp_tbi)
            .map { meta, file, index -> tuple([vartype: meta.vartype] + [tag: "TP_comp"] + [id: "rtgtools"], file, index) }
            .set { vcf_tp_comp }

        tagged_variants = tagged_variants.mix(vcf_fn)
        tagged_variants = tagged_variants.mix(vcf_fp)
        tagged_variants = tagged_variants.mix(vcf_tp_base)
        tagged_variants = tagged_variants.mix(vcf_tp_comp)
    }

    if (params.method.contains('happy')){

        HAPPY_HAPPY(
            input_ch.map { it -> tuple(it[0], it[1], it[3], it[5], []) },
            fasta,
            fai,
            [[],[]],
            [[],[]],
            [[],[]]
        )
        versions = versions.mix(HAPPY_HAPPY.out.versions)

        HAPPY_HAPPY.out.summary_csv
            .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "happy"], file) }
            .groupTuple()
            .set{ report}
        summary_reports = summary_reports.mix(report)
    }
    emit:
    versions
    summary_reports
    tagged_variants

}
