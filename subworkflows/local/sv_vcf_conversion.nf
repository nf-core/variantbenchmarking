import groovy.io.FileType

//
// SV_VCF_CONVERSIONS: SUBWORKFLOW TO apply tool spesific conversions
//

params.options = [:]

include { SVYNC                   } from '../../modules/nf-core/svync'
include { BGZIP_TABIX             } from '../../modules/local/bgzip_tabix'
include { VARIANT_EXTRACTOR       } from '../../modules/local/variant_extractor'
include { BCFTOOLS_SORT           } from '../../modules/nf-core/bcftools/sort'

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


    emit:
    vcf_ch
    versions
}
