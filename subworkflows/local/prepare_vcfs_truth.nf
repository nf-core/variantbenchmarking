//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//


include { BCFTOOLS_NORM              } from '../../modules/nf-core/bcftools/norm'
include { TABIX_TABIX                } from '../../modules/nf-core/tabix/tabix'
include { VCF_REHEADER_SAMPLENAME    } from '../local/vcf_reheader_samplename'
include { VCF_VARIANT_DEDUPLICATION  } from '../local/vcf_variant_deduplication'
include { LIFTOVER_VCFS_TRUTH        } from '../local/liftover_vcfs_truth'


workflow PREPARE_VCFS_TRUTH {
    take:
    truth_ch        // channel: [val(meta), vcf]
    high_conf_ch    // channel: [val(meta), bed]
    fasta           // reference channel [val(meta), ref.fa]
    fai             // reference channel [val(meta), ref.fa.fai]
    chain           // reference channel [val(meta), chain.gz]
    liftover_genome // reference channel [val(meta), ref.fa]
    rename_chr      // reference channel [val(meta), chrlist.txt]

    main:

    versions = Channel.empty()
    bed_high_conf = Channel.empty()

    // if liftover option is set convert truth files
    if (params.liftover){

        LIFTOVER_VCFS_TRUTH(
            truth_ch,
            high_conf_ch,
            liftover_genome,
            chain,
            rename_chr
        )
        versions = versions.mix(LIFTOVER_VCFS_TRUTH.out.versions.first())
        truth_ch = LIFTOVER_VCFS_TRUTH.out.vcf_ch
        bed_high_conf = LIFTOVER_VCFS_TRUTH.out.bed_ch
    }

    // Reheader sample name for truth file - using meta.caller
    VCF_REHEADER_SAMPLENAME(
        truth_ch,
        fai
    )
    versions = versions.mix(VCF_REHEADER_SAMPLENAME.out.versions.first())
    vcf_ch   = VCF_REHEADER_SAMPLENAME.out.ch_vcf

    if (params.preprocess.contains("normalization")){

        // multi-allelic variants will be splitted.
        BCFTOOLS_NORM(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_NORM.out.versions.first())

        // index vcf file
        TABIX_TABIX(
            BCFTOOLS_NORM.out.vcf
        )
        versions = versions.mix(TABIX_TABIX.out.versions)

        BCFTOOLS_NORM.out.vcf.join(TABIX_TABIX.out.tbi, by:0)
                            .set{vcf_ch}
    }
    if (params.preprocess.contains("deduplication")){

        // Deduplicates variants at the same position test
        VCF_VARIANT_DEDUPLICATION(
            vcf_ch,
            fasta
        )
        vcf_ch = VCF_VARIANT_DEDUPLICATION.out.ch_vcf
        versions = versions.mix(VCF_VARIANT_DEDUPLICATION.out.versions)
    }

    emit:
    vcf_ch
    bed_high_conf
    versions
}
