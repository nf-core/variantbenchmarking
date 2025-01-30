//
// CNV_GERMLINE_BENCHMARK: SUBWORKFLOW FOR CNV GERMLINE VARIANTS
//

include { WITTYER                          } from '../../../modules/nf-core/wittyer'
include { TABIX_BGZIP as TABIX_BGZIP_QUERY } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_BGZIP as TABIX_BGZIP_TRUTH } from '../../../modules/nf-core/tabix/bgzip'

workflow CNV_GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),test_vcf,test_index,truth_vcf,truth_index, bed]

    main:

    versions        = Channel.empty()
    summary_reports = Channel.empty()

    // CNV benchmarking is only possible with wittyer now!

    TABIX_BGZIP_QUERY(
        input_ch.map{ meta, vcf, _tbi, _truth_vcf, _truth_tbi, _bed ->
            [ meta, vcf ]
        }
    )
    versions = versions.mix(TABIX_BGZIP_QUERY.out.versions.first())

    TABIX_BGZIP_TRUTH(
        input_ch.map{ meta, _vcf, _tbi, truth_vcf, _truth_tbi, _bed ->
            [ meta, truth_vcf ]
        }
    )
    versions = versions.mix(TABIX_BGZIP_TRUTH.out.versions.first())

    input_ch.map{ meta, _vcf, _tbi, _truth_vcf, _truth_tbi, bedfile ->
            [ meta, bedfile ]
        }
        .set { bed }

    // RUN WITTYER for benchmarking
    WITTYER(
        TABIX_BGZIP_QUERY.out.output
            .join(TABIX_BGZIP_TRUTH.out.output, failOnMismatch: true, failOnDuplicate: true)
            .join(bed, failOnMismatch: true, failOnDuplicate: true)
    )
    versions = versions.mix(WITTYER.out.versions.first())

    WITTYER.out.report
        .map { meta, file ->
            tuple([vartype: meta.vartype] + [benchmark_tool: "wittyer"], file)
        }
        .groupTuple()
        .set{ report }
    summary_reports = summary_reports.mix(report)


    emit:
    summary_reports  // channel: [val(meta), summary]
    versions         // channel: [versions.yml]
}
