//
// SV_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SV GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_PHAB           } from '../../modules/local/truvari_phab'                  addParams( options: params.options )
include { TRUVARI_BENCH          } from '../../modules/nf-core/truvari/bench'          addParams( options: params.options )
include { SVANALYZER_SVBENCHMARK } from '../../modules/nf-core/svanalyzer/svbenchmark' addParams( options: params.options )

workflow SV_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),test_vcf,test_index,truth_vcf,truth_index, bed]
    fasta     // reference channel [val(meta), ref.fa]
    fai       // reference channel [val(meta), ref.fa.fai]

    main:

    versions=Channel.empty()
    summary_reports=Channel.empty()
    tagged_variants=Channel.empty()

    // SV benchmarking

    if (params.method.contains('truvari')){

        if(params.sv_standardization.contains('harmonize')){
            //
            // TRUVARI: TRUVARI_PHAB
            //
            TRUVARI_PHAB(
                input_ch,
                fasta,
                fai
            )
        }
        //
        // MODULE: TRUVARI_BENCH
        //
        TRUVARI_BENCH(
            input_ch,
            fasta,
            fai
        )
        versions = versions.mix(TRUVARI_BENCH.out.versions)

        TRUVARI_BENCH.out.summary
            .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "truvari"], file) }
            .groupTuple()
            .set { report }

        summary_reports = summary_reports.mix(report)

        TRUVARI_BENCH.out.fn_vcf
            .map { meta, file -> tuple([vartype: meta.vartype] + [tag: "FN"] + [id: "truvari"], file) }
            .set { vcf_fn }

        TRUVARI_BENCH.out.fp_vcf
            .map { meta, file -> tuple([vartype: meta.vartype] + [tag: "FP"] + [id: "truvari"], file) }
            .set { vcf_fp }

        TRUVARI_BENCH.out.tp_base_vcf
            .map { meta, file -> tuple([vartype: meta.vartype] + [tag: "TP_base"] + [id: "truvari"], file) }
            .set { vcf_tp_base }

        TRUVARI_BENCH.out.tp_comp_vcf
            .map { meta, file -> tuple([vartype: meta.vartype] + [tag: "TP_comp"] + [id: "truvari"], file) }
            .set { vcf_tp_comp }

        tagged_variants = tagged_variants.mix(vcf_fn)
        tagged_variants = tagged_variants.mix(vcf_fp)
        tagged_variants = tagged_variants.mix(vcf_tp_base)
        tagged_variants = tagged_variants.mix(vcf_tp_comp)

    }

    if (params.method.contains('svanalyzer')){
        //
        // MODULE: SVANALYZER_SVBENCHMARK
        //
        // slower than truvari
        SVANALYZER_SVBENCHMARK(
            input_ch,
            fasta,
            fai
            )
        versions = versions.mix(SVANALYZER_SVBENCHMARK.out.versions)

        SVANALYZER_SVBENCHMARK.out.report
            .map { meta, file -> tuple([vartype: meta.vartype] + [benchmark_tool: "svbenchmark"], file) }
            .groupTuple()
            .set{ report}

        summary_reports = summary_reports.mix(report)

        SVANALYZER_SVBENCHMARK.out.fns
            .map { meta, file -> tuple([vartype: meta.vartype] + [tag: "FN"] + [id: "svbenchmark"], file) }
            .set { vcf_fn }

        SVANALYZER_SVBENCHMARK.out.fps
            .map { meta, file -> tuple([vartype: meta.vartype] + [tag: "FP"] + [id: "svbenchmark"], file) }
            .set { vcf_fp }
        tagged_variants = tagged_variants.mix(vcf_fn)
        tagged_variants = tagged_variants.mix(vcf_fp)

    }

    emit:
    tagged_variants
    summary_reports
    versions
}
