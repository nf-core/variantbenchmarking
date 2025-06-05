//
// WITTYER_BENCHMARK: SUBWORKFLOW FOR BENCHMARKING SV VARIANTS WITH WITTYER
//

include { WITTYER    } from '../../../modules/nf-core/wittyer'
include { TABIX_BGZIP as TABIX_BGZIP_QUERY  } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_BGZIP as TABIX_BGZIP_TRUTH  } from '../../../modules/nf-core/tabix/bgzip'

workflow WITTYER_BENCHMARK {
    take:
    input_ch  // channel: [val(meta), test_vcf, test_index, truth_vcf, truth_index, regionsbed, targets_bed ]

    main:

    versions        = Channel.empty()

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

    emit:
    report         // channel: [val(meta), reports]
    versions     // channel: [val(meta), versions.yml]
}
