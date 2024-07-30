//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

include { BCFTOOLS_NORM     } from '../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_REHEADER } from '../../modules/nf-core/bcftools/reheader'
include { VCF_VARIANT_DEDUPLICATION      } from '../local/vcf_variant_deduplication'
include { TABIX_TABIX as TABIX_TABIX_1   } from '../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_TABIX_2   } from '../../modules/nf-core/tabix/tabix'


workflow PREPARE_VCFS_TRUTH {
    take:
    truth_ch    // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions=Channel.empty()

    //
    // BCFTOOLS_REHEADER
    //
    BCFTOOLS_REHEADER(
        truth_ch.map{meta, file -> tuple( meta, file, [], [])},
        fai
        )
    versions = versions.mix(BCFTOOLS_REHEADER.out.versions)

    //
    // TABIX_TABIX
    //
    TABIX_TABIX_1(
        BCFTOOLS_REHEADER.out.vcf
    )
    versions = versions.mix(TABIX_TABIX_1.out.versions)
    BCFTOOLS_REHEADER.out.vcf.join(TABIX_TABIX_1.out.tbi, by:0)
                            .set{vcf_ch}


    if (params.preprocess.contains("normalization")){
        //
        // MODULE:  BCFTOOLS_NORM
        //
        // Normalize test
        // multi-allelic variants will be splitted.
        BCFTOOLS_NORM(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_NORM.out.versions)
        //
        // TABIX_BGZIPTABIX
        //
        // index vcf file
        TABIX_TABIX_2(
            BCFTOOLS_NORM.out.vcf
        )
        versions = versions.mix(TABIX_TABIX_2.out.versions)

        BCFTOOLS_NORM.out.vcf.join(TABIX_TABIX_2.out.tbi, by:0)
                            .set{vcf_ch}
    }
    if (params.preprocess.contains("deduplication")){
        //
        // VCF_VARIANT_DEDUPLICATION
        //
        // Deduplicates variants at the same position test
        VCF_VARIANT_DEDUPLICATION(
            vcf_ch,
            fasta
        )
        vcf_ch = VCF_VARIANT_DEDUPLICATION.out.ch_vcf
        }

    emit:
    vcf_ch
    versions
}
