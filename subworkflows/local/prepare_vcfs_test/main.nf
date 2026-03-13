//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

include { VCF_VARIANT_DEDUPLICATION    } from '../../local/vcf_variant_deduplication'
include { VCF_VARIANT_FILTERING        } from '../../local/vcf_variant_filtering'
include { SPLIT_SMALL_VARIANTS_TEST    } from '../../local/split_small_variants_test'
include { LIFTOVER_VCFS                } from '../../local/liftover_vcfs'
include { BCFTOOLS_NORM                } from '../../../modules/nf-core/bcftools/norm'
include { FIX_VCF_PREFIX               } from '../../../modules/local/custom/fix_vcf_prefix'
include { RTGTOOLS_SVDECOMPOSE         } from '../../../modules/nf-core/rtgtools/svdecompose'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CONTIGS      } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_NORM as BCFTOOLS_SPLIT_MULTI       } from '../../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_QUERY} from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FILTERMISSING} from '../../../modules/nf-core/bcftools/view'
include { GAWK as ADD_GT_STRELKA                      } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIPTABIX as TABIX_BGZIPTABIX_GT     } from '../../../modules/nf-core/tabix/bgziptabix'

workflow PREPARE_VCFS_TEST {
    take:
    test_ch     // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]
    chain       // reference channel [val(meta), chain.gz]
    rename_chr  // reference channel [val(meta), chrlist.txt]
    dictionary  // reference channel [val(meta), genome.dict]

    main:

    versions = channel.empty()

    // branch out test samples with metadata liftover is true
    test_ch.branch{input ->
        def meta = input[0]
        liftover: meta.liftover
        other: true}.set{vcf}

    vcf_ch = channel.empty()

    if (params.liftover.contains("test")){

        // apply liftover test vcfs
        LIFTOVER_VCFS(
            vcf.liftover,
            channel.empty(),
            fasta,
            chain,
            rename_chr,
            dictionary
        )
        vcf_ch = vcf_ch.mix(LIFTOVER_VCFS.out.vcf_ch)
    }
    vcf_ch = vcf_ch.mix(vcf.other)

    // if prefix of chromosomes needs to be fixed
    vcf_ch.branch{ input ->
        def meta = input[0]
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

    // rename sample name
    BCFTOOLS_REHEADER_QUERY(

        vcf_ch.map{ meta, file ->
            [ meta, file, [], [] ]
        },
        fai
    )

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

        BCFTOOLS_VIEW_CONTIGS.out.vcf.join(BCFTOOLS_VIEW_CONTIGS.out.tbi, by:0)
                            .set{vcf_ch}
    }
    if (params.preprocess.contains("split_multiallelic")){

        // Split -any- multi-allelic variants
        BCFTOOLS_SPLIT_MULTI(
            vcf_ch,
            fasta
        )

        BCFTOOLS_SPLIT_MULTI.out.vcf.join(BCFTOOLS_SPLIT_MULTI.out.tbi, by:0)
                            .set{vcf_ch}
    }

    if (params.include_expression != null || params.exclude_expression != null || params.min_sv_size > 0 || params.max_sv_size != -1 || params.min_allele_freq != -1 || params.min_num_reads != -1 ){

        // Filters variants and SVs with given parameters
        VCF_VARIANT_FILTERING(
            vcf_ch
        )
        vcf_ch = VCF_VARIANT_FILTERING.out.vcf_ch
    }

    if (params.preprocess.contains("deduplicate")){
        // Deduplicate variants at the same position test
        VCF_VARIANT_DEDUPLICATION(
            vcf_ch,
            fasta
        )
        vcf_ch = VCF_VARIANT_DEDUPLICATION.out.ch_vcf
    }

    if (params.preprocess.contains("normalize")){
        // Turn on left alignment and normalization
        BCFTOOLS_NORM(
            vcf_ch,
            fasta
        )
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, by:0)
                            .set{vcf_ch}
    }

    if (params.analysis.contains("somatic")){
        // somatic specific preparations
        if (params.variant_type == "small"){
            SPLIT_SMALL_VARIANTS_TEST(
                vcf_ch
            )
            vcf_ch = SPLIT_SMALL_VARIANTS_TEST.out.out_vcf_ch
        }
    }

    if (params.sv_standardization.contains("svdecompose")){
        RTGTOOLS_SVDECOMPOSE(
            vcf_ch
        )
        vcf_ch = RTGTOOLS_SVDECOMPOSE.out.vcf.join(RTGTOOLS_SVDECOMPOSE.out.index)
    }

    // filters out ./. or 0/0 or non-somatic genotypes
    BCFTOOLS_VIEW_FILTERMISSING(
        vcf_ch,
        [],
        [],
        []
    )

     // branch out test samples with missing GT field (e.g. strelka) to add GT field using gawk
    ch_branched_vcf = BCFTOOLS_VIEW_FILTERMISSING.out.vcf.join(BCFTOOLS_VIEW_FILTERMISSING.out.tbi)
        .branch { meta, vcf, tbi ->
            def is_rtg = params.method?.contains("rtgtools")
            def is_strelka_manta = ['strelka', 'manta'].contains(meta.caller.toLowerCase())
            def is_somatic = params.analysis == "somatic"

            needs_gt: is_strelka_manta && is_somatic && is_rtg
            ok:       true
        }

    // Add GT field using
    ADD_GT_STRELKA(
        ch_branched_vcf.needs_gt.map{ meta, vcf, _tbi -> tuple(meta, vcf) },
        [],
        false
    )

    TABIX_BGZIPTABIX_GT(
        ADD_GT_STRELKA.out.output
    )

    vcf_ch = TABIX_BGZIPTABIX_GT.out.gz_index.mix(ch_branched_vcf.ok)

    emit:
    vcf_ch   // channel: [val(meta), vcf.gz, tbi]
    versions // channel: [versions.yml]
}
