//
// SV_VCF_CONVERSIONS: SUBWORKFLOW to apply tool spesific conversions
//

include { SVYNC                   } from '../../../modules/nf-core/svync'
include { BGZIP_TABIX             } from '../../../modules/local/bgzip/tabix'
include { VARIANT_EXTRACTOR       } from '../../../modules/local/custom/variant_extractor'
include { SVTK_STANDARDIZE        } from '../../../modules/nf-core/svtk/standardize/main.nf'
include { BCFTOOLS_SORT as BCFTOOLS_SORT1 } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_SORT as BCFTOOLS_SORT2 } from '../../../modules/nf-core/bcftools/sort'

workflow SV_VCF_CONVERSIONS {
    take:
    input_ch    // channel: [val(meta), vcf]
    fasta       // reference channel [val(meta), ref.fa]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:
    versions   = Channel.empty()

    if (params.sv_standardization.contains("variant_extractor")){
        // uses VariantExtractor to homogenize variants
        VARIANT_EXTRACTOR(
            input_ch,
            fasta,
            fai
        )
        versions = versions.mix(VARIANT_EXTRACTOR.out.versions.first())

        // sort vcf
        BCFTOOLS_SORT1(
            VARIANT_EXTRACTOR.out.output
        )
        versions = versions.mix(BCFTOOLS_SORT1.out.versions.first())
        input_ch = BCFTOOLS_SORT1.out.vcf

    }

    if (params.sv_standardization.contains("svtk_standardize")){

        out_vcf_ch = Channel.empty()

        supported_callers2 = ["delly", "melt", "manta", "wham"]

        input_ch
            .branch{ meta, vcf->
                def caller = meta.caller
                def supported = supported_callers2.contains(caller)
                if(!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported. Skipping standardization...")
                }
                tool:  supported
                other: !supported
            }
            .set{input}

        SVTK_STANDARDIZE(
            input.tool,
            fai.map{meta, file -> file}
        )
        versions = versions.mix(SVTK_STANDARDIZE.out.versions.first())

        BCFTOOLS_SORT2(
            SVTK_STANDARDIZE.out.standardized_vcf
        )
        versions = versions.mix(BCFTOOLS_SORT2.out.versions.first())

        out_vcf_ch.mix(
                BCFTOOLS_SORT2.out.vcf,
                input.other
            ).set{input_ch}

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
        supported_callers = ["delly", "dragen", "gridss", "manta", "smoove"]

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
