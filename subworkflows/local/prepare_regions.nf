//
// PREPARE_STRATIFICATIONS: SUBWORKFLOW TO PREPARE BED FILES - HIGH CONFIDENCE AND OTHER LEVELS
//

params.options = [:]

include { MAIN_CHROMS         } from '../../modules/local/main_chroms.nf'            addParams( options: params.options )
include { EXTRACT_MAIN        } from '../../modules/local/extract_main.nf'           addParams( options: params.options )


workflow PREPARE_REGIONS {
    take:
    ref       // reference channel [ref.fa, ref.fa.fai]
    high_conf

    main:

    versions=Channel.empty()

    ref.map { it -> tuple([id: it[0].baseName], it[1]) }
            .set{fasta}

    // this is not working!
        
    // get contig file including only main chroms
    MAIN_CHROMS(
        fasta
    )
    main_chroms = MAIN_CHROMS.out.sizes
    versions = versions.mix(MAIN_CHROMS.out.versions)

    high_conf.map { it -> tuple([id: it[0].baseName], it[0]) }
            .set{bed}
    bed.view()
    EXTRACT_MAIN(
        bed
    )
    chr_list = EXTRACT_MAIN.out.chr_list

    emit:
    main_chroms
    chr_list
    versions
}
