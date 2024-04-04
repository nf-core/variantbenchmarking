//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

params.options = [:]

include { BGZIP_TABIX      } from '../../modules/local/bgzip_tabix.nf'       addParams( options: params.options )
include { BCFTOOLS_VIEW    } from '../../modules/local/bcftools_view'        addParams( options: params.options )
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix'   addParams( options: params.options )
include { BCFTOOLS_NORM    } from '../../modules/nf-core/bcftools/norm'      addParams( options: params.options )
include { TABIX_TABIX      } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { VCF_VARIANT_DEDUPLICATION                    } from '../local/vcf_variant_deduplication'      addParams( options: params.options )
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_TRUTH } from '../../modules/nf-core/bcftools/reheader' addParams( options: params.options )

workflow PREPARE_VCFS_TRUTH {
    take:
    truth_ch    // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    versions=Channel.empty()

    //
    // PREPARE_VCFS
    //

    // BGZIP if needed and index truth
    BGZIP_TABIX(
        truth_ch
    )
    versions = versions.mix(BGZIP_TABIX.out.versions)
    vcf_ch = BGZIP_TABIX.out.gz_tbi

    // Reheader needed to standardize sample names
    ch_samples = Channel.of(["samples.txt", params.sample,"truth"])
                    .collectFile(newLine:false)

    vcf_ch.combine(ch_samples)
            .map{it -> tuple( it[0], it[1],[],it[3])}
            .set{input_ch}
    //
    // BCFTOOLS_REHEADER
    //
    // Add "truth" to test sample
    BCFTOOLS_REHEADER_TRUTH(
        input_ch,
        fai
        )
    versions = versions.mix(BCFTOOLS_REHEADER_TRUTH.out.versions)

    //
    // TABIX_BGZIPTABIX
    //
    // bgzip and index vcf file
    TABIX_BGZIPTABIX(
        BCFTOOLS_REHEADER_TRUTH.out.vcf
    )
    vcf_ch = TABIX_BGZIPTABIX.out.gz_tbi

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
        // ndex vcf file
        TABIX_TABIX(
            BCFTOOLS_NORM.out.vcf
        )
        versions = versions.mix(TABIX_TABIX.out.versions)

        BCFTOOLS_NORM.out.vcf.join(TABIX_TABIX.out.tbi, by:0)
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
