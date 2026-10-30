//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

include { VCF_VARIANT_DEDUPLICATION  } from '../../local/vcf_variant_deduplication'
include { LIFTOVER_VCFS              } from '../../local/liftover_vcfs'
include { BCFTOOLS_NORM              } from '../../../modules/nf-core/bcftools/norm'
include { RTGTOOLS_SVDECOMPOSE       } from '../../../modules/nf-core/rtgtools/svdecompose'
include { BCFTOOLS_NORM as BCFTOOLS_SPLIT_MULTI        } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_TRUTH } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTERMISSING } from '../../../modules/nf-core/bcftools/view'

workflow PREPARE_VCFS_TRUTH {
    take:
    truth_ch        // channel: [val(meta), vcf]
    regions_bed_ch  // channel: [bed]
    targets_bed_ch  // channel: [bed]
    fasta           // reference channel [val(meta), ref.fa]
    fai             // reference channel [val(meta), ref.fa.fai]
    chain           // reference channel [val(meta), chain.gz]
    rename_chr      // reference channel [val(meta), chrlist.txt]
    dictionary      // reference channel [val(meta), genome.dict]

    main:

    // if liftover option is set convert truth files
    if (params.liftover.contains("truth")){
        LIFTOVER_VCFS(
            truth_ch,
            regions_bed_ch,
            targets_bed_ch,
            fasta,
            chain,
            rename_chr,
            dictionary
        )
        truth_ch = LIFTOVER_VCFS.out.vcf_ch
        regions_bed_ch = LIFTOVER_VCFS.out.bed_ch.map{ _meta, bed -> [bed]}
        targets_bed_ch = LIFTOVER_VCFS.out.targets_ch.map{ _meta, bed -> [bed]}

    }

    // rename sample name
    BCFTOOLS_REHEADER_TRUTH(
        truth_ch.map{ meta, vcf ->
            [ meta, vcf, [], [] ]
        },
        fai
    )

    BCFTOOLS_REHEADER_TRUTH.out.vcf.join(BCFTOOLS_REHEADER_TRUTH.out.index)
        .set{vcf_ch}

    if (params.preprocess.contains("split_multiallelic")){
        // Split -any- multi-allelic variants
        BCFTOOLS_SPLIT_MULTI(
            vcf_ch,
            fasta
        )
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
    }

    if (params.preprocess.contains("normalize")){

        // Turn on left alignment and m\normalization
        BCFTOOLS_NORM(
            vcf_ch,
            fasta
        )
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, by:0)
                            .set{vcf_ch}
    }

    if (params.sv_standardization.contains("svdecompose")){
        RTGTOOLS_SVDECOMPOSE(
            vcf_ch
        )
        vcf_ch = RTGTOOLS_SVDECOMPOSE.out.vcf.join(RTGTOOLS_SVDECOMPOSE.out.index)
    }

    if (!(params.enable_missing_genotypes?.contains("truth"))) {
        // filters out ./. or 0/0 or non-somatic genotypes
        BCFTOOLS_VIEW_FILTERMISSING(
            vcf_ch,
            [],
            [],
            []
        )
        vcf_ch = BCFTOOLS_VIEW_FILTERMISSING.out.vcf.join(BCFTOOLS_VIEW_FILTERMISSING.out.tbi)
    }

    emit:
    vcf_ch         // channel: [val(meta), vcf, tbi]
    regions_bed_ch // channel: [val(meta), bed]
    targets_bed_ch // channel: [val(meta), bed]
}
