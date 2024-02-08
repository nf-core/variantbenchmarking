//
// VCF_CONVERSIONS: SUBWORKFLOW TO apply tool spesific conversions
//

params.options = [:]

include { MANTA_CONVERTINVERSION  } from '../../modules/nf-core/manta/convertinversion'  addParams( options: params.options )
include { GRIDSS_ANNOTATION       } from '../../modules/local/gridss_annotation'         addParams( options: params.options )
include { SVYNC                   } from '../../modules/nf-core/svync'                   addParams( options: params.options )
include { AWK_SORT                } from '../../modules/local/awk_sort.nf'               addParams( options: params.options )

workflow VCF_CONVERSIONS {
    take:
    input_ch    // channel: [val(meta),val(meta2), vcf, config.yml]
    ref         // reference channel [ref.fa, ref.fa.fai]

    main:

    versions=Channel.empty()

    //
    // MODULE: AWK_SORT
    //
    // sort and index input test files

    AWK_SORT(
        input_ch.map{it -> tuple(it[0], it[1], it[2])}
    ) 
    versions = versions.mix(AWK_SORT.out.versions)
    vcf_ch = AWK_SORT.out.vcf

    //
    // MODULE: SVYNC
    //
    // 
    if(params.standardization){ 

        input_ch.map{it -> tuple(it[0], it[1], it[3])}
            .combine(vcf_ch, by:1)
            .map{it -> tuple(it[1], it[0], it[4], it[5], it[2])}
            .set{snd_ch} 
        
        SVYNC(
            snd_ch
        )
        vcf_ch = SVYNC.out.vcf
    }
    

    // Check tool spesific conversions
    if(params.bnd_to_inv){
        //
        // MANTA_CONVERTINVERSION
        //
        //NOTE: should also work for dragen
        MANTA_CONVERTINVERSION(
            vcf_ch,
            ref
        )
        versions = versions.mix(MANTA_CONVERTINVERSION.out.versions)
        vcf_ch = MANTA_CONVERTINVERSION.out.vcf_tabi

       // https://github.com/srbehera/DRAGEN_Analysis/blob/main/convertInversion.py

    }


    if (params.gridss_annotate){
        //
        // GRIDSS_ANNOTATION
        //
        // https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
        // GRIDSS simple event annotation
        GRIDSS_ANNOTATION(
            vcf_ch,
            ref
        )
        versions = versions.mix(GRIDSS_ANNOTATION.out.versions)
        vcf_ch = GRIDSS_ANNOTATION.out.vcf
    }

   // https://github.com/EUCANCan/variant-extractor/blob/main/examples/vcf_to_csv.py

    emit:
    vcf_ch
    versions
}
