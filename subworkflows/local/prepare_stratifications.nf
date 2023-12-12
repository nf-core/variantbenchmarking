//
// PREPARE_STRATIFICATIONS: SUBWORKFLOW TO PREPARE BED FILES - HIGH CONFIDENCE AND OTHER LEVELS
//

params.options = [:]

include { BEDTOOLS_INTERSECT  } from '../../modules/nf-core/bedtools/intersect'      addParams( options: params.options )
include { MAIN_CHROMS         } from '../../modules/local/main_chroms.nf'            addParams( options: params.options )


workflow PREPARE_STRATIFICATIONS {
    take:
    ref       // reference channel [ref.fa, ref.fa.fai]

    main:

    versions=Channel.empty()

    ref.map { it -> tuple([id: it[0].baseName], it[1]) }
            .set{fasta}

    // get contig file including only main chroms
    MAIN_CHROMS(
        fasta
    )
    main_chroms = MAIN_CHROMS.out.sizes
    versions = versions.mix(MAIN_CHROMS.out.versions)

    emit:
    main_chroms
    versions
}
