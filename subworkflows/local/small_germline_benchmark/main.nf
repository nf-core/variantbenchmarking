//
// SMALL_GERMLINE_BENCHMARK: SUBWORKFLOW FOR SMALL GERMLINE VARIANTS
//

include { RTGTOOLS_BENCHMARK  } from '../../../subworkflows/local/rtgtools_benchmark'
include { HAPPY_BENCHMARK     } from '../../../subworkflows/local/happy_benchmark'


workflow SMALL_GERMLINE_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targetsbed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    sdf                // reference channel [val(meta), sdf]
    falsepositive_bed  // reference channel [val(meta), bed]
    stratification_bed // reference channel [val(meta), bed files]
    stratification_tsv // reference channel [val(meta), tsv]

    main:

    versions        = Channel.empty()
    summary_reports = Channel.empty()
    tagged_variants = Channel.empty()

    if (params.method.contains('rtgtools')){

        RTGTOOLS_BENCHMARK(
            input_ch,
            fasta,
            fai,
            sdf
        )
        versions        = versions.mix(RTGTOOLS_BENCHMARK.out.versions.first())
        summary_reports = summary_reports.mix(RTGTOOLS_BENCHMARK.out.summary_reports)
        tagged_variants = tagged_variants.mix(RTGTOOLS_BENCHMARK.out.tagged_variants)

    }

    if (params.method.contains('happy')){

        HAPPY_BENCHMARK(
            input_ch,
            fasta,
            fai,
            falsepositive_bed,
            stratification_bed,
            stratification_tsv

        )
        versions        = versions.mix(HAPPY_BENCHMARK.out.versions)
        summary_reports = summary_reports.mix(HAPPY_BENCHMARK.out.summary_reports)
        tagged_variants = tagged_variants.mix(HAPPY_BENCHMARK.out.tagged_variants)

    }
    emit:
    summary_reports // channel: [val(meta), reports]
    tagged_variants // channel: [val(meta), vcfs]
    versions        // channel: [versions.yml]

}
