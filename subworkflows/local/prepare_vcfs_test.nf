//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

params.options = [:]

include { BCFTOOLS_REHEADER   } from '../../modules/nf-core/bcftools/reheader'  addParams( options: params.options )
include { BCFTOOLS_RENAME_CHR } from '../../modules/local/bcftools_rename_chr'  addParams( options: params.options )
include { BCFTOOLS_VIEW       } from '../../modules/local/bcftools_view'        addParams( options: params.options )
include { BCFTOOLS_NORM as BCFTOOLS_NORM_1  } from '../../modules/nf-core/bcftools/norm'  addParams( options: params.options )
include { BCFTOOLS_NORM as BCFTOOLS_NORM_2  } from '../../modules/nf-core/bcftools/norm'  addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_1    } from '../../modules/nf-core/tabix/tabix'    addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_2    } from '../../modules/nf-core/tabix/tabix'    addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_3    } from '../../modules/nf-core/tabix/tabix'    addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_4    } from '../../modules/nf-core/tabix/tabix'    addParams( options: params.options )
include { BGZIP_TABIX as BGZIP_TABIX_1      } from '../../modules/local/bgzip_tabix'      addParams( options: params.options )
include { BGZIP_TABIX as BGZIP_TABIX_2      } from '../../modules/local/bgzip_tabix'      addParams( options: params.options )

workflow PREPARE_VCFS_TEST {
    take:
    input_ch    // channel: [val(meta),val(meta2), vcf]
    ref         // reference channel [ref.fa, ref.fa.fai]
    rename_chr  // channel path(rename chromosomes)
    main_chroms // channel path(chrom sizes)
    chr_list

    main:

    versions=Channel.empty()

    //ref.map { it -> tuple([id: it[0].baseName], it[1]) }
    //    .set{fasta}
    //
    // PREPARE_VCFS
    //
    // Reheader needed to standardize sample names
    BCFTOOLS_REHEADER(
        input_ch,
        ref,
        params.sample
        )
    versions = versions.mix(BCFTOOLS_REHEADER.out.versions)


    BGZIP_TABIX_1(
        BCFTOOLS_REHEADER.out.vcf
    )
    reheader_ch = BGZIP_TABIX_1.out.gz_tbi

    //TABIX_TABIX_1(
    //    BCFTOOLS_REHEADER.out.vcf
    //)
    //versions = versions.mix(TABIX_TABIX_1.out.versions)

    //BCFTOOLS_REHEADER.out.vcf.join(TABIX_TABIX_1.out.tbi, by:1)
    //                    .map{it -> tuple(it[1], it[0], it[2], it[4])}
    //                    .set{reheader_ch}
    reheader_ch.view()
    // 1 -> chr1 or chr1 -> 1
    BCFTOOLS_RENAME_CHR(
        reheader_ch,
        rename_chr
    )
    versions = versions.mix(BCFTOOLS_RENAME_CHR.out.versions)

    TABIX_TABIX_2(
        BCFTOOLS_RENAME_CHR.out.vcf
    )
    versions = versions.mix(TABIX_TABIX_2.out.versions)

    BCFTOOLS_RENAME_CHR.out.vcf.join(TABIX_TABIX_2.out.tbi, by:1)
                        .map{it -> tuple(it[1], it[0], it[2], it[4])}
                        .set{vcf_ch}
    //
    // BCFTOOLS_VIEW
    //
    // To filter out contigs!
    BCFTOOLS_VIEW(
        vcf_ch
    )
    versions = versions.mix(BCFTOOLS_VIEW.out.versions)

    BGZIP_TABIX_2(
        BCFTOOLS_VIEW.out.vcf
    )
    versions = versions.mix(BGZIP_TABIX_2.out.versions)
    vcf_ch   = BGZIP_TABIX_2.out.gz_tbi

    if (params.preprocess.contains("normalization")){
        //
        // BCFTOOLS_NORM
        //
        // Breaks down -any- multi-allelic variants 
        BCFTOOLS_NORM_1(
            vcf_ch,
            ref,
            [[],[]]
        )
        versions = versions.mix(BCFTOOLS_NORM_1.out.versions)

        TABIX_TABIX_3(
            BCFTOOLS_NORM_1.out.vcf
        )
        versions = versions.mix(TABIX_TABIX_3.out.versions)

        BCFTOOLS_NORM_1.out.vcf.join(TABIX_TABIX_3.out.tbi, by:1)
                            .map{it -> tuple( it[1], it[0], it[2], it[4])}
                            .set{vcf_ch}
    }

    if (params.preprocess.contains("deduplication")){
        //
        // BCFTOOLS_NORM
        //
        // Deduplicates variants at the same position test
        BCFTOOLS_NORM_2(
            vcf_ch,
            ref,
            [[],[]]
        )
        versions = versions.mix(BCFTOOLS_NORM_2.out.versions)

        TABIX_TABIX_4(
            BCFTOOLS_NORM_2.out.vcf
        )
        versions = versions.mix(TABIX_TABIX_4.out.versions)

        BCFTOOLS_NORM_2.out.vcf.join(TABIX_TABIX_4.out.tbi, by:1)
                            .map{it -> tuple( it[1], it[0], it[2], it[4])}
                            .set{vcf_ch}
    }
    emit:
    vcf_ch
    versions
}
