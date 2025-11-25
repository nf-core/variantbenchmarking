//
// PREPARE_REFERENCES: SUBWORKFLOW TO PREPARE REFERENCES AND DICTINARIES
//

include { PICARD_CREATESEQUENCEDICTIONARY   } from '../../../modules/nf-core/picard/createsequencedictionary'
include { RTGTOOLS_FORMAT                   } from '../../../modules/nf-core/rtgtools/format/main'

workflow PREPARE_REFERENCES {
    take:
    fasta_ch   // reference channel [val(meta), genome.fasta]
    dictionary // reference channel [val(meta), genome.dict]
    sdf        // reference channel [val(meta), genome.sdf]

    main:

    versions = Channel.empty()

    //prepare dict file for liftover of vcf files
    if (!params.dictionary && (params.method.contains("concordance") || params.liftover ) ){
        PICARD_CREATESEQUENCEDICTIONARY(
            fasta_ch
        )
        dictionary = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
        versions = versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)
    }

    if (!params.sdf && params.method.contains("rtgtools")){

        // Use rtgtools format to generate sdf file if necessary
        RTGTOOLS_FORMAT(
            fasta_ch.map { meta, file -> [ meta, file, [], [] ] }
        )
        versions = versions.mix(RTGTOOLS_FORMAT.out.versions)
        sdf = RTGTOOLS_FORMAT.out.sdf
    }

    emit:
    dictionary   // reference channel [val(meta), genome.dict]
    sdf          // reference channel [val(meta), genome.sdf]
    versions     // channel: [val(meta), versions.yml]

}
