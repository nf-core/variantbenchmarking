import groovy.io.FileType

//
// SV_VCF_CONVERSIONS: SUBWORKFLOW TO apply tool spesific conversions
//

params.options = [:]

include { MANTA_CONVERTINVERSION  } from '../../modules/nf-core/manta/convertinversion'  addParams( options: params.options )
include { GRIDSS_ANNOTATION       } from '../../modules/local/gridss_annotation'         addParams( options: params.options )
include { SVYNC                   } from '../../modules/nf-core/svync'                   addParams( options: params.options )
include { BGZIP_TABIX             } from '../../modules/local/bgzip_tabix'               addParams( options: params.options )
include { VARIANT_EXTRACTOR       } from '../../modules/local/variant_extractor'         addParams( options: params.options )
include { BCFTOOLS_SORT           } from '../../modules/nf-core/bcftools/sort'           addParams( options: params.options )

workflow SV_VCF_CONVERSIONS {
    take:
    input_ch    // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]
    svync_yaml  // reference channel [yamls]

    main:
    versions   = Channel.empty()

    if (params.sv_standardization.contains("homogenize")){
        // uses VariantExtractor to homogenize variants

        VARIANT_EXTRACTOR(
            input_ch,
            fasta,
            fai
        )
        versions = versions.mix(VARIANT_EXTRACTOR.out.versions)

        // sort vcf
        BCFTOOLS_SORT(
            VARIANT_EXTRACTOR.out.output
        )
        versions = versions.mix(BCFTOOLS_SORT.out.versions)
        input_ch = BCFTOOLS_SORT.out.vcf

    }

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
        supported_callers = []
        new File("${projectDir}/assets/svync").eachFileRecurse (FileType.FILES) { supported_callers << it.baseName.replace(".yaml", "") }

        svync_yaml
            .map { yamls ->
                [yamls.collect { it.baseName }, yamls.collectEntries { [it.baseName, it] }]
            }
            .set { custom_svync_configs }

        vcf_ch
            .combine(custom_svync_configs)
            .branch{ meta, vcf, tbi, custom_callers, custom_configs ->
                def caller = meta.id
                def supported = custom_callers.contains(caller) || supported_callers.contains(caller)
                if(!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported. Skipping standardization...")
                }
                tool:  supported
                    return [ meta, vcf, tbi, custom_configs[caller] ?: [] ]
                other: !supported
                    return [ meta, vcf, tbi ]
            }
            .set{input}

        input.tool
            .map { meta, vcf, tbi, custom_config ->
                config = custom_config ?: file("${projectDir}/assets/svync/${meta.id}.yaml", checkIfExists:true)
                [ meta, vcf, tbi, config ]
            }
            .set {svync_ch}

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
    // https://github.com/EUCANCan/variant-extractor/blob/main/examples/vcf_to_csv.py'


    emit:
    vcf_ch
    versions
}
