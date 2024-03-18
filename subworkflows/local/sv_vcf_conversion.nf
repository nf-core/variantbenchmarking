//
// SV_VCF_CONVERSIONS: SUBWORKFLOW TO apply tool spesific conversions
//

params.options = [:]

include { MANTA_CONVERTINVERSION  } from '../../modules/nf-core/manta/convertinversion'  addParams( options: params.options )
include { GRIDSS_ANNOTATION       } from '../../modules/local/gridss_annotation'         addParams( options: params.options )
include { SVYNC                   } from '../../modules/nf-core/svync'                   addParams( options: params.options )
include { BGZIP_TABIX             } from '../../modules/local/bgzip_tabix'               addParams( options: params.options )

workflow SV_VCF_CONVERSIONS {
    take:
    input_ch    // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]
    svync_yaml  // yaml configs

    main:
    versions   = Channel.empty()

    //
    // MODULE: BGZIP_TABIX
    //
    // zip and index input test files

    BGZIP_TABIX(
        input_ch
    )
    versions = versions.mix(BGZIP_TABIX.out.versions)
    vcf_ch = BGZIP_TABIX.out.gz_tbi

    //
    // MODULE: SVYNC
    //
    //
    if(params.sv_standardization.contains("standardization")){
        out_vcf_ch = Channel.empty()

        vcf_ch.branch{
            tool:  it[0].id == "delly" || it[0].id == "gridss" || it[0].id == "manta" || it[0].id == "smoove"
            other: true}
            .set{input}

        svync_yaml.flatten()
            .map { it -> tuple([id: it.getSimpleName(), vartype: "sv"], it) }
            .set{config}

        input.tool
            .combine(config, by:0)
            .set{svync_ch}

        SVYNC(
            svync_ch
        )
        out_vcf_ch = out_vcf_ch.mix(SVYNC.out.vcf)
        out_vcf_ch = out_vcf_ch.mix(input.other)
        vcf_ch     = out_vcf_ch.map{it -> tuple(it[0], it[1], it[2])}
    }

    // Check tool spesific conversions
    if(params.sv_standardization.contains("bnd_to_inv")){
        out_vcf_ch = Channel.empty()

        vcf_ch.branch{
            tool:  it[0].id == "manta" || it[0].id == "dragen"
            other: true}
            .set{input}
        //
        // MANTA_CONVERTINVERSION
        //
        // NOTE: should also work for dragen
        // Not working now!!!!!

        MANTA_CONVERTINVERSION(
            input.tool.map{it -> tuple(it[0], it[1])},
            fasta
        )
        versions = versions.mix(MANTA_CONVERTINVERSION.out.versions)

        out_ch = MANTA_CONVERTINVERSION.out.vcf.join(MANTA_CONVERTINVERSION.out.tbi)
        out_vcf_ch = out_vcf_ch.mix(out_ch)
        out_vcf_ch = out_vcf_ch.mix(input.other)
        vcf_ch     = out_vcf_ch

        // https://github.com/srbehera/DRAGEN_Analysis/blob/main/convertInversion.py

    }

    if (params.sv_standardization.contains("gridss_annotate")){
        out_vcf_ch = Channel.empty()

        vcf_ch.branch{
            tool:  it[0].id == "gridss"
            other: true}
            .set{input}

        //
        // GRIDSS_ANNOTATION
        //
        // https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
        // GRIDSS simple event annotation
        GRIDSS_ANNOTATION(
            input.tool,
            fasta,
            fai
        )
        versions = versions.mix(GRIDSS_ANNOTATION.out.versions)

        out_vcf_ch = out_vcf_ch.mix(GRIDSS_ANNOTATION.out.vcf)
        out_vcf_ch = out_vcf_ch.mix(input.other)
        vcf_ch     = out_vcf_ch
    }
    // https://github.com/EUCANCan/variant-extractor/blob/main/examples/vcf_to_csv.py

    emit:
    vcf_ch
    versions
}
