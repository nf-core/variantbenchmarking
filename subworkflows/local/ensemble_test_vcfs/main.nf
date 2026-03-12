
//
// ENSEMBLE_TEST_VCFS: SUBWORKFLOW to ENSEMBLE TEST VCFS and PREPARE TRUTH VCF.
//

include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_SMALL } from '../../../modules/nf-core/tabix/bgziptabix'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_GT    } from '../../../modules/nf-core/tabix/bgziptabix'
include { BCFTOOLS_MERGE as BCFTOOLS_ENSEMBLE } from '../../../modules/nf-core/bcftools/merge'
include { SURVIVOR_MERGE as SURVIVOR_ENSEMBLE } from '../../../modules/nf-core/survivor/merge'
include { BCFTOOLS_VIEW as FILTER_MAJORITY    } from '../../../modules/nf-core/bcftools/view'
include { GAWK as REFORMAT_TRUTH              } from '../../../modules/nf-core/gawk'
include { GAWK as REFORMAT_TRUTH_SV           } from '../../../modules/nf-core/gawk'
include { GAWK as INJECT_MISSING_GT           } from '../../../modules/nf-core/gawk'
include { BCFTOOLS_SORT as BCFTOOLS_SORT_SV   } from '../../../modules/nf-core/bcftools/sort'
include { TABIX_BGZIP as TABIX_BGZIP_UNZIP    } from '../../../modules/nf-core/tabix/bgzip'

workflow ENSEMBLE_TEST_VCFS {
    take:
    test_vcfs        // channel: [val(meta), vcf.gz, index]
    fasta           // reference channel [val(meta), ref.fa]
    fai             // reference channel [val(meta), ref.fa.fai]

    main:

    // if the benchmarking method is rtgtools, missing GT field is already filled in VCF preperation step, so no need to inject missing GT field
    test_vcfs.branch { meta, vcf, index ->
        def is_rtg = params.method?.contains("rtgtools")
        def is_strelka_manta = ['strelka', 'manta'].contains(meta.caller.toLowerCase())
        def is_somatic = params.analysis == "somatic"

        missing_gt: is_strelka_manta && is_somatic && !is_rtg
        other:       true
    }.set { branched_vcfs }

    // Run the injection process ONLY on  strelka
    INJECT_MISSING_GT(
        branched_vcfs.missing_gt.map { meta, vcf, index -> tuple(meta, vcf) },
        [],
        false
    )

    if (params.variant_type == "small" | params.variant_type == "snv" | params.variant_type == "indel"){

        TABIX_BGZIPTABIX_GT(
            INJECT_MISSING_GT.out.output
        )

        ch_ready_for_merge = branched_vcfs.other.mix(
            TABIX_BGZIPTABIX_GT.out.gz_index
        )

        // drop meta information from vcf test samples
        ch_test_vcfs = ch_ready_for_merge.map { meta, vcf, index ->
            [ [id: 'truth'], vcf, index ]
        }

        // merge small variants
        BCFTOOLS_ENSEMBLE(
            ch_test_vcfs.groupTuple(),
            fasta,
            fai,
            [[],[]]
        )

        FILTER_MAJORITY(
            BCFTOOLS_ENSEMBLE.out.vcf.join(BCFTOOLS_ENSEMBLE.out.index),
            [], //regions
            [],
            [],
        )

        REFORMAT_TRUTH(
            FILTER_MAJORITY.out.vcf,
            [],
            false
        )

        TABIX_BGZIPTABIX_SMALL(
            REFORMAT_TRUTH.out.output
        )

        truth_vcf = TABIX_BGZIPTABIX_SMALL.out.gz_index

    }
    else{
        // Merge SV variants

        TABIX_BGZIP_UNZIP(
            branched_vcfs.other.map { meta, vcf, _index -> tuple(meta, vcf) }
        )

        TABIX_BGZIP_UNZIP.out.output.mix(INJECT_MISSING_GT.out.output)
            .set{vcf_ch}

        // drop meta information from vcf test samples
        ch_test_vcfs = vcf_ch.map { meta, vcf ->
            [ [id: 'truth'], vcf ]
        }

        // Dynamically calculate the minimum callers based on the number of files
        // Merge Benchmark SVs from different tools
        SURVIVOR_ENSEMBLE(
            ch_test_vcfs.groupTuple(),
            10000,
            1,
            params.ensemble_truth ,
            0,
            0,
            params.min_sv_size
        )

        REFORMAT_TRUTH_SV(
            SURVIVOR_ENSEMBLE.out.vcf,
            [],
            false
        )

        BCFTOOLS_SORT_SV(
            REFORMAT_TRUTH_SV.out.output
        )
        truth_vcf = BCFTOOLS_SORT_SV.out.vcf.join(BCFTOOLS_SORT_SV.out.tbi)

    }
    emit:
    truth_vcf    // channel: [val(meta), vcf.gz, index]

}
