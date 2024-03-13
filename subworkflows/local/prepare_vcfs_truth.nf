//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

params.options = [:]

include { BGZIP_TABIX      } from '../../modules/local/bgzip_tabix.nf'       addParams( options: params.options )
include { BCFTOOLS_VIEW    } from '../../modules/local/bcftools_view'        addParams( options: params.options )
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix'   addParams( options: params.options )
include { BCFTOOLS_NORM as BCFTOOLS_NORM_1 } from '../../modules/nf-core/bcftools/norm'      addParams( options: params.options )
include { BCFTOOLS_NORM as BCFTOOLS_NORM_2 } from '../../modules/nf-core/bcftools/norm'      addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_1   } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_2   } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_3   } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_TRUTH } from '../../modules/nf-core/bcftools/reheader'  addParams( options: params.options )

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
    ch_samples = Channel.of(["samples.txt", params.sample,"_truth"])
                    .collectFile(newLine:false)

    vcf_ch.combine(ch_samples)
            .map{it -> tuple( it[0], it[1],[],it[3])}
            .set{input_ch}

    BCFTOOLS_REHEADER_TRUTH(
        input_ch,
        fai
        )
    versions = versions.mix(BCFTOOLS_REHEADER_TRUTH.out.versions)

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
        BCFTOOLS_NORM_1(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_NORM_1.out.versions)

        TABIX_TABIX_1(
            BCFTOOLS_NORM_1.out.vcf
        )
        versions = versions.mix(TABIX_TABIX_1.out.versions)

        BCFTOOLS_NORM_1.out.vcf.join(TABIX_TABIX_1.out.tbi, by:0)
                            .set{vcf_ch}
    }
    if (params.preprocess.contains("deduplication")){
        //
        // MODULE:  BCFTOOLS_NORM
        //
        // Deduplicate variants at the same position
        BCFTOOLS_NORM_2(
            vcf_ch,
            fasta
        )
        versions = versions.mix(BCFTOOLS_NORM_2.out.versions)

        TABIX_TABIX_2(
            BCFTOOLS_NORM_2.out.vcf
        )
        versions = versions.mix(TABIX_TABIX_2.out.versions)

        BCFTOOLS_NORM_2.out.vcf.join(TABIX_TABIX_2.out.tbi, by:0)
                            .set{vcf_ch}
        }


    emit:
    vcf_ch
    versions
}
