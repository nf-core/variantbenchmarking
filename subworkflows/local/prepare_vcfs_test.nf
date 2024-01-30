//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

params.options = [:]

include { BCFTOOLS_REHEADER   } from '../../modules/nf-core/bcftools/reheader'  addParams( options: params.options )
include { BCFTOOLS_NORM       } from '../../modules/nf-core/bcftools/norm'      addParams( options: params.options )
include { BCFTOOLS_SORT       } from '../../modules/nf-core/bcftools/sort'      addParams( options: params.options )
include { AWK_SORT            } from '../../modules/local/awk_sort'             addParams( options: params.options )
include { BGZIP_TABIX         } from '../../modules/local/bgzip_tabix'          addParams( options: params.options )
include { BCFTOOLS_RENAME_CHR } from '../../modules/local/bcftools_rename_chr'  addParams( options: params.options )
include { TABIX_TABIX as TABIX_TABIX_1  } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { TABIX_TABIX as TABIX_TABIX_2  } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { TABIX_TABIX as TABIX_TABIX_3  } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )

workflow PREPARE_VCFS_TEST {
    take:
    input_ch    // channel: [val(meta), vcf]
    ref         // reference channel [ref.fa, ref.fa.fai]
    rename_chr  // channel path(rename chromosomes)
    main_chroms // channel path(chrom sizes)

    main:

    versions=Channel.empty()

    ref.map { it -> tuple([id: it[0].baseName], it[1]) }
        .set{fasta}
    //
    // PREPARE_VCFS
    //
    // Reheader needed to standardize sample names
    BCFTOOLS_REHEADER(
        input_ch,
        fasta
    )
    versions = versions.mix(BCFTOOLS_REHEADER.out.versions)

    AWK_SORT(
        BCFTOOLS_REHEADER.out.vcf
    )
    //BCFTOOLS_SORT.out.vcf.join(BCFTOOLS_SORT.out.index)
    //                    .map{it -> tuple(it[0], it[1], it[2], it[4])}  
    //                    .set{reheader_ch}    

    BGZIP_TABIX(
        AWK_SORT.out.vcf
    )
    reheader_ch = BGZIP_TABIX.out.gz_tbi

    //TABIX_TABIX_1(
   //     BCFTOOLS_REHEADER.out.vcf
    //)
    //versions = versions.mix(TABIX_TABIX_1.out.versions)

    //BCFTOOLS_REHEADER.out.vcf.join(TABIX_TABIX_1.out.tbi)
    //                    .map{it -> tuple(it[0], it[1], it[2], it[4])}
    //                    .set{reheader_ch}

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
                        .set{renamed_ch}
    // Normalize test
    BCFTOOLS_NORM(
        renamed_ch,
        ref,
        [[],[]]
    )
    versions = versions.mix(BCFTOOLS_NORM.out.versions)

    TABIX_TABIX_3(
        BCFTOOLS_NORM.out.vcf
    )
    versions = versions.mix(TABIX_TABIX_3.out.versions)

    BCFTOOLS_NORM.out.vcf.join(TABIX_TABIX_3.out.tbi, by:1)
                        .map{it -> tuple( it[1], it[0], it[2], it[4])}
                        .set{vcf_ch}

    emit:
    vcf_ch
    versions
}
