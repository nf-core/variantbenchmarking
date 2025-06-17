//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//


include { VCF_VARIANT_DEDUPLICATION  } from '../../local/vcf_variant_deduplication'
include { LIFTOVER_VCFS              } from '../../local/liftover_vcfs'
include { BCFTOOLS_NORM              } from '../../../modules/nf-core/bcftools/norm'
include { PUBLISH_PROCESSED_VCF      } from '../../../modules/local/custom/publish_processed_vcf'
include { BCFTOOLS_NORM as BCFTOOLS_SPLIT_MULTI       } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_TRUTH} from '../../../modules/local/bcftools/reheader'
include { RTGTOOLS_SVDECOMPOSE    } from '../../../modules/nf-core/rtgtools/svdecompose'


workflow PREPARE_VCFS_TRUTH {
    take:
    truth_ch        // channel: [val(meta), vcf]
    high_conf_ch    // channel: [bed]
    fasta           // reference channel [val(meta), ref.fa]
    fai             // reference channel [val(meta), ref.fa.fai]
    chain           // reference channel [val(meta), chain.gz]
    rename_chr      // reference channel [val(meta), chrlist.txt]
    dictionary      // reference channel [val(meta), genome.dict]

    main:

    versions = Channel.empty()

    // if liftover option is set convert truth files
    if (params.liftover.contains("truth")){

        LIFTOVER_VCFS(
            truth_ch,
            high_conf_ch,
            fasta,
            chain,
            rename_chr,
            dictionary
        )
        versions = versions.mix(LIFTOVER_VCFS.out.versions)
        truth_ch = LIFTOVER_VCFS.out.vcf_ch
        high_conf_ch = LIFTOVER_VCFS.out.bed_ch.map{ _meta, bed -> [bed]}
    }

    // Reheader sample name for truth file - using meta.caller
    // rename sample name
    BCFTOOLS_REHEADER_TRUTH(
        truth_ch.map{ meta, file ->
            [ meta, file, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_TRUTH.out.versions)

    BCFTOOLS_REHEADER_TRUTH.out.vcf.join(BCFTOOLS_REHEADER_TRUTH.out.index)
        .set{vcf_ch}

    if (params.preprocess.contains("split_multiallelic")){

        // Split -any- multi-allelic variants
        BCFTOOLS_SPLIT_MULTI(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_SPLIT_MULTI.out.versions)

        BCFTOOLS_SPLIT_MULTI.out.vcf.join(BCFTOOLS_SPLIT_MULTI.out.tbi, by:0)
                            .set{vcf_ch}
    }

    if (params.preprocess.contains("deduplicate")){

        // Deduplicates variants at the same position test
        VCF_VARIANT_DEDUPLICATION(
            vcf_ch,
            fasta
        )
        vcf_ch = VCF_VARIANT_DEDUPLICATION.out.ch_vcf
        versions = versions.mix(VCF_VARIANT_DEDUPLICATION.out.versions)
    }

    if (params.preprocess.contains("normalize")){

        // Turn on left alignment and m\normalization
        BCFTOOLS_NORM(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_NORM.out.versions)

        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, by:0)
                            .set{vcf_ch}
    }

    if (params.sv_standardization.contains("svdecompose")){
        RTGTOOLS_SVDECOMPOSE(
            vcf_ch
        )
        versions = versions.mix(RTGTOOLS_SVDECOMPOSE.out.versions)
        vcf_ch = RTGTOOLS_SVDECOMPOSE.out.vcf.join(RTGTOOLS_SVDECOMPOSE.out.index)
    }

    PUBLISH_PROCESSED_VCF(
        vcf_ch
    )
    versions = versions.mix(PUBLISH_PROCESSED_VCF.out.versions)

    emit:
    vcf_ch       // channel: [val(meta), vcf, tbi]
    high_conf_ch // channel: [val(meta), bed]
    versions     // channel: [versions.yml]
}
