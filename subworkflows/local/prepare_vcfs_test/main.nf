//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

include { VCF_VARIANT_DEDUPLICATION    } from '../../local/vcf_variant_deduplication'
include { VCF_VARIANT_FILTERING        } from '../../local/vcf_variant_filtering'
include { SPLIT_SMALL_VARIANTS_TEST    } from '../../local/split_small_variants_test'
include { LIFTOVER_VCFS                } from '../../local/liftover_vcfs'
include { BCFTOOLS_NORM                } from '../../../modules/nf-core/bcftools/norm'
include { FIX_VCF_PREFIX               } from '../../../modules/local/custom/fix_vcf_prefix'
include { PUBLISH_PROCESSED_VCF        } from '../../../modules/local/custom/publish_processed_vcf'
include { RTGTOOLS_SVDECOMPOSE         } from '../../../modules/nf-core/rtgtools/svdecompose'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CONTIGS      } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_NORM as BCFTOOLS_SPLIT_MULTI       } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_QUERY} from '../../../modules/local/bcftools/reheader'


workflow PREPARE_VCFS_TEST {
    take:
    test_ch     // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]
    chain       // reference channel [val(meta), chain.gz]
    rename_chr  // reference channel [val(meta), chrlist.txt]
    dictionary  // reference channel [val(meta), genome.dict]

    main:

    versions = Channel.empty()

    // branch out test samples with metadata liftover is true
    test_ch.branch{
        def meta = it[0]
        liftover: meta.liftover
        other: true}.set{vcf}

    vcf_ch = Channel.empty()

    if (params.liftover.contains("test")){

        // apply liftover test vcfs
        LIFTOVER_VCFS(
            vcf.liftover,
            Channel.empty(),
            fasta,
            chain,
            rename_chr,
            dictionary
        )
        versions = versions.mix(LIFTOVER_VCFS.out.versions.first())
        vcf_ch = vcf_ch.mix(LIFTOVER_VCFS.out.vcf_ch)
    }
    vcf_ch = vcf_ch.mix(vcf.other)

    // if prefix of chromosomes needs to be fixed
    vcf_ch.branch{
        def meta = it[0]
        prefix: meta.fix_prefix
        other: true}.set{fix}

    vcf_ch = Channel.empty()

    // fix vcf chromosome prefix according to reference genome
    FIX_VCF_PREFIX(
        fix.prefix,
        rename_chr
    )
    versions = versions.mix(FIX_VCF_PREFIX.out.versions.first())
    vcf_ch = vcf_ch.mix(FIX_VCF_PREFIX.out.vcf,fix.other)


    // Add "query" to test sample
    // rename sample name
    BCFTOOLS_REHEADER_QUERY(

        vcf_ch.map{ meta, file ->
            [ meta, file, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER_QUERY.out.versions.first())

    BCFTOOLS_REHEADER_QUERY.out.vcf.join(BCFTOOLS_REHEADER_QUERY.out.index)
        .set{vcf_ch}

    if (params.preprocess.contains("filter_contigs")){
        // filter out extra contigs!
        BCFTOOLS_VIEW_CONTIGS(
            vcf_ch,
            [],
            [],
            []
        )
        versions = versions.mix(BCFTOOLS_VIEW_CONTIGS.out.versions.first())

        BCFTOOLS_VIEW_CONTIGS.out.vcf.join(BCFTOOLS_VIEW_CONTIGS.out.tbi, by:0)
                            .set{vcf_ch}
    }
    if (params.preprocess.contains("split_multiallelic")){

        // Split -any- multi-allelic variants
        BCFTOOLS_SPLIT_MULTI(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_SPLIT_MULTI.out.versions.first())

        BCFTOOLS_SPLIT_MULTI.out.vcf.join(BCFTOOLS_SPLIT_MULTI.out.tbi, by:0)
                            .set{vcf_ch}
    }

    if (params.include_expression != null || params.exclude_expression != null || params.min_sv_size > 0 || params.max_sv_size != -1 || params.min_allele_freq != -1 || params.min_num_reads != -1 ){

        // Filters variants and SVs with given parameters
        VCF_VARIANT_FILTERING(
            vcf_ch
        )
        vcf_ch = VCF_VARIANT_FILTERING.out.vcf_ch
        versions = versions.mix(VCF_VARIANT_FILTERING.out.versions.first())
    }

    if (params.preprocess.contains("deduplicate")){

        // Deduplicate variants at the same position test
        VCF_VARIANT_DEDUPLICATION(
            vcf_ch,
            fasta
        )
        vcf_ch = VCF_VARIANT_DEDUPLICATION.out.ch_vcf
        versions = versions.mix(VCF_VARIANT_DEDUPLICATION.out.versions.first())

    }

    if (params.preprocess.contains("normalize")){

        // Turn on left alignment and normalization
        BCFTOOLS_NORM(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_NORM.out.versions.first())

        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, by:0)
                            .set{vcf_ch}
    }

    if (params.analysis.contains("somatic")){

        // somatic specific preparations
        if (params.variant_type == "small"){
            SPLIT_SMALL_VARIANTS_TEST(
                vcf_ch
            )
            versions = versions.mix(SPLIT_SMALL_VARIANTS_TEST.out.versions.first())
            vcf_ch = SPLIT_SMALL_VARIANTS_TEST.out.out_vcf_ch
        }

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

    emit:
    vcf_ch   // channel: [val(meta), vcf.gz, tbi]
    versions // channel: [versions.yml]
}
