//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

include { VCF_VARIANT_DEDUPLICATION              } from '../local/vcf_variant_deduplication'
include { VCF_VARIANT_FILTERING                  } from '../local/vcf_variant_filtering'
include { SUBSAMPLE_SOMATIC_VCFS_TEST            } from '../local/subsample_somatic_vcfs_test'
include { SPLIT_SMALL_VARIANTS_TEST              } from '../local/split_small_variants_test'
include { BGZIP_TABIX                            } from '../../modules/local/bgzip_tabix'
include { BCFTOOLS_NORM                          } from '../../modules/nf-core/bcftools/norm'
include { BCFTOOLS_REHEADER                      } from '../../modules/nf-core/bcftools/reheader'
include { TABIX_BGZIPTABIX                       } from '../../modules/nf-core/tabix/bgziptabix'
include { TABIX_TABIX as TABIX_TABIX_1           } from '../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_TABIX_2           } from '../../modules/nf-core/tabix/tabix'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_CONTIGS } from '../../modules/nf-core/bcftools/view'


workflow PREPARE_VCFS_TEST {
    take:
    input_ch    // channel: [val(meta),vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions = Channel.empty()
    //
    // MODULE: BGZIP_TABIX
    //
    // zip and index input test files
    BGZIP_TABIX(
        input_ch
    )
    versions = versions.mix(BGZIP_TABIX.out.versions.first())
    vcf_ch = BGZIP_TABIX.out.gz_tbi

    if (params.analysis.contains("somatic")){

        out_vcf_ch = Channel.empty()
        // subsample multisample vcf if necessary

        vcf_ch.branch{
                def meta = it[0]
                sample: meta.subsample != null
                other:  true
            }
            .set{vcf}

        SUBSAMPLE_SOMATIC_VCFS_TEST(
            vcf.sample
        )
        versions = versions.mix(SUBSAMPLE_SOMATIC_VCFS_TEST.out.versions)
        out_vcf_ch = out_vcf_ch.mix(SUBSAMPLE_SOMATIC_VCFS_TEST.out.vcf_ch)
        out_vcf_ch = out_vcf_ch.mix(vcf.other)
        vcf_ch = out_vcf_ch
    }

    //
    // BCFTOOLS_REHEADER
    //
    // Add "query" to test sample
    BCFTOOLS_REHEADER(
        vcf_ch.map{ meta, vcf, tbi -> 
            [ meta, vcf, [], [] ]
        },
        fai
    )
    versions = versions.mix(BCFTOOLS_REHEADER.out.versions.first())

    TABIX_TABIX_1(
        BCFTOOLS_REHEADER.out.vcf
    )
    versions = versions.mix(TABIX_TABIX_1.out.versions.first())
    BCFTOOLS_REHEADER.out.vcf
        .join(TABIX_TABIX_1.out.tbi, failOnDuplicate: true, failOnMismatch:true)
        .set{vcf_ch}

    if (params.preprocess.contains("filter_contigs")){
        //
        // BCFTOOLS_VIEW
        //
        // To filter out contigs!
        BCFTOOLS_VIEW_CONTIGS(
            vcf_ch,
            [],
            [],
            []
        )
        versions = versions.mix(BCFTOOLS_VIEW_CONTIGS.out.versions.first())

        TABIX_BGZIPTABIX(
            BCFTOOLS_VIEW_CONTIGS.out.vcf
        )
        versions = versions.mix(TABIX_BGZIPTABIX.out.versions.first())
        vcf_ch   = TABIX_BGZIPTABIX.out.gz_tbi
    }
    if (params.preprocess.contains("normalization")){
        //
        // BCFTOOLS_NORM
        //
        // Breaks down -any- multi-allelic variants
        BCFTOOLS_NORM(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_NORM.out.versions.first())

        TABIX_TABIX_2(
            BCFTOOLS_NORM.out.vcf
        )
        versions = versions.mix(TABIX_TABIX_2.out.versions.first())
        BCFTOOLS_NORM.out.vcf
            .join(TABIX_TABIX_2.out.tbi, failOnDuplicate: true, failOnMismatch: true)
            .set{vcf_ch}
    }

    if (params.include_expression != null | params.exclude_expression != null | params.min_sv_size > 0 | params.max_sv_size != -1 | params.min_allele_freq != -1 | params.min_num_reads != -1 ){

        //
        // SUBWORKFLOW: VCF_VARIANT_FILTERING
        //
        // Filters SVs with given paramaters
        VCF_VARIANT_FILTERING(
            vcf_ch
        )
        vcf_ch = VCF_VARIANT_FILTERING.out.vcf_ch
        versions = versions.mix(VCF_VARIANT_FILTERING.out.versions)
    }

    if (params.preprocess.contains("deduplication")){
        //
        // SUBWORKFLOW: VCF_VARIANT_DEDUPLICATION
        //
        // Deduplicates variants at the same position test
        VCF_VARIANT_DEDUPLICATION(
            vcf_ch,
            fasta
        )
        vcf_ch = VCF_VARIANT_DEDUPLICATION.out.ch_vcf
        versions = versions.mix(VCF_VARIANT_DEDUPLICATION.out.versions)

    }

    if (params.analysis.contains("somatic")){

        // somatic spesific preperations
        vcf_ch.branch{
                def meta = it[0]
                small: meta.vartype == "small"
                other: true
            }
            .set{vcf}

        out_vcf_ch = Channel.empty()

        SPLIT_SMALL_VARIANTS_TEST(
            vcf.small
        )
        versions = versions.mix(SPLIT_SMALL_VARIANTS_TEST.out.versions)
        out_vcf_ch = out_vcf_ch.mix(
            SPLIT_SMALL_VARIANTS_TEST.out.out_vcf_ch,
            vcf.other
        )
        vcf_ch = out_vcf_ch
    }

    emit:
    vcf_ch
    versions
}
