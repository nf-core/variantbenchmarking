//
// SV_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SV GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_PHAB           } from '../../modules/local/truvari_phab'                  addParams( options: params.options )
include { TRUVARI_BENCH          } from '../../modules/nf-core/truvari/bench'          addParams( options: params.options )
include { SVANALYZER_SVBENCHMARK } from '../../modules/nf-core/svanalyzer/svbenchmark' addParams( options: params.options )
include { WITTYER                } from '../../modules/nf-core/wittyer'                  addParams( options: params.options )
include { VCFDIST                } from '../../modules/local/vcfdist'                  addParams( options: params.options )

workflow SV_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),test_vcf,test_index,truth_vcf,truth_index, bed]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions=Channel.empty()
    summary_reports=Channel.empty()

    // SV benchmarking

    if (params.method.contains('truvari')){

        if(params.harmonize){
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

    }

    if (params.method.contains('wittyer')){

        //
        // MODULE: WITTYER
        //
        // BIG Advantage: reports by variant type
        // Able to report CNV
        WITTYER(
            input_ch
        )
        versions = versions.mix(WITTYER.out.versions)
    }

    if (params.method.contains('vcfdist')){
        //
        // MODULE: VCFDIST
        //
        VCFDIST(
            input_ch,
            fasta,
            fai
        )
        versions = versions.mix(VCFDIST.out.versions)
    }

    emit:
    summary_reports
    versions
}
