//
// SOMATIC: SUBWORKFLOW FOR SMALL SOMATIC VARIANTS
//

include { SOMPY_BENCHMARK     } from '../../../subworkflows/local/sompy_benchmark'
include { RTGTOOLS_BENCHMARK  } from '../../../subworkflows/local/rtgtools_benchmark'


workflow SMALL_SOMATIC_BENCHMARK {
    take:
    input_ch           // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index,  regionsbed, targets_bed ]
    fasta              // reference channel [val(meta), ref.fa]
    fai                // reference channel [val(meta), ref.fa.fai]
    sdf                // reference channel [val(meta), sdf]
    falsepositive_bed  // reference channel [val(meta), bed]
    ambiguous_beds     // reference channel [val(meta), bed]

    main:

    versions            = Channel.empty()
    summary_reports     = Channel.empty()
    tagged_variants     = Channel.empty()
    tagged_variants_csv = Channel.empty()

    if (params.method.contains('sompy')){
        // apply sompy for small somatic variant benchmarking
        SOMPY_BENCHMARK(
            input_ch,
            fasta,
            fai,
            falsepositive_bed,
            ambiguous_beds
        )
        versions        = versions.mix(SOMPY_BENCHMARK.out.versions)
        summary_reports = summary_reports.mix(SOMPY_BENCHMARK.out.summary_reports)
        tagged_variants = tagged_variants_csv.mix(SOMPY_BENCHMARK.out.tagged_variants_csv)
    }

    if (params.method.contains('rtgtools')){

        RTGTOOLS_BENCHMARK(
            input_ch,
            fasta,
            fai,
            sdf
        )
        versions        = versions.mix(RTGTOOLS_BENCHMARK.out.versions)
        summary_reports = summary_reports.mix(RTGTOOLS_BENCHMARK.out.summary_reports)
        tagged_variants = tagged_variants.mix(RTGTOOLS_BENCHMARK.out.tagged_variants)

    }

    emit:
    summary_reports     // channel: [val(meta), reports]
    tagged_variants     // channel: [val(meta), vcfs]
    tagged_variants_csv // channel: [val(meta), csvs]
    versions            // channel: [versions.yml]
}
