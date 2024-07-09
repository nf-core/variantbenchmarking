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
    if(params.sv_standardization.contains("svync")){
        out_vcf_ch = Channel.empty()
        supported_callers = ["delly", "dragen", "gridss", "manta", "delly", "smoove"]

        vcf_ch
            .branch{ meta, vcf, tbi ->
                def caller = meta.id
                def supported = supported_callers.contains(caller)
                if(!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported. Skipping standardization...")
                }
                tool:  supported
                    return [ meta, vcf, tbi]
                other: !supported
                    return [ meta, vcf, tbi ]
            }
            .set{input}


        input.tool
            .map { meta, vcf, tbi ->
                [ meta, vcf, tbi, file("${projectDir}/assets/svync/${meta.id}.yaml", checkIfExists:true) ]
            }
            .set {svync_ch}

        svync_ch.view()

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
