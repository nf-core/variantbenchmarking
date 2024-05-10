//
// REPORT_BENCHMARK_STATISTICS: BENCHMARK REPORTS
//

params.options = [:]

include { BGZIP_TABIX      } from '../../modules/local/bgzip_tabix.nf'       addParams( options: params.options )
include { BCFTOOLS_VIEW    } from '../../modules/local/bcftools_view'      addParams( options: params.options )

workflow REPORT_BENCHMARK_STATISTICS {
    take:
    truth_ch    // channel: [val(meta), vcf]
    ref         // reference channel [ref.fa, ref.fa.fai]
    main_chroms // channel: path(chrom sizes)

    main:

    versions=Channel.empty()

// Check tool spesific conversions

    // https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
    // GRIDSS simple event annotation

    emit:
    vcf_ch
    versions
}
