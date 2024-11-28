import groovy.io.FileType

//
// SV_VCF_CONVERSIONS: SUBWORKFLOW to apply tool spesific conversions
//

include { SVYNC                   } from '../../modules/nf-core/svync'
include { BGZIP_TABIX             } from '../../modules/local/bgzip/tabix'
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
        versions = versions.mix(VARIANT_EXTRACTOR.out.versions.first())

        // sort vcf
        BCFTOOLS_SORT(
            VARIANT_EXTRACTOR.out.output
        )
        versions = versions.mix(BCFTOOLS_SORT.out.versions.first())
        input_ch = BCFTOOLS_SORT.out.vcf

    }

    // zip and index input test files
    BGZIP_TABIX(
        input_ch
    )
    versions = versions.mix(BGZIP_TABIX.out.versions.first())
    vcf_ch = BGZIP_TABIX.out.gz_tbi

    // RUN SVYNC tool to reformat SV callers
    if(params.sv_standardization.contains("svync")){
        out_vcf_ch = Channel.empty()
        supported_callers = ["delly", "dragen", "gridss", "manta", "delly", "smoove"]

        vcf_ch
            .branch{ meta, vcf, tbi ->
                def caller = meta.caller
                def supported = supported_callers.contains(caller)
                if(!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported. Skipping standardization...")
                }
                tool:  supported
                    return [ meta, vcf, tbi]
                other: !supported
                    return [ meta, vcf ]
            }
            .set{input}


        input.tool
            .map { meta, vcf, tbi ->
                [ meta, vcf, tbi, file("${projectDir}/assets/svync/${meta.caller}.yaml", checkIfExists:true) ]
            }
            .set {svync_ch}

        SVYNC(
            svync_ch
        )
        versions = versions.mix(SVYNC.out.versions.first())
        out_vcf_ch.mix(
                SVYNC.out.vcf,
                input.other
            )
            .map{
                def meta = it[0]
                def vcf = it[1]
                [ meta, vcf ]
            }
            .set { vcf_ch }
    }

    emit:
    vcf_ch   // channel: [val(meta), vcf]
    versions // channel: [versions.yml]
}
