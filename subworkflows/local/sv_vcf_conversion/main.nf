//
// SV_VCF_CONVERSIONS: SUBWORKFLOW to apply tool spesific conversions
//

include { SVYNC                   } from '../../../modules/nf-core/svync'
include { TABIX_BGZIPTABIX        } from '../../../modules/nf-core/tabix/bgziptabix'
include { VARIANTEXTRACTOR        } from '../../../modules/nf-core/variantextractor'
include { SVTK_STANDARDIZE        } from '../../../modules/nf-core/svtk/standardize'
include { RTGTOOLS_SVDECOMPOSE    } from '../../../modules/nf-core/rtgtools/svdecompose'
include { BCFTOOLS_SORT as BCFTOOLS_SORT1 } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_SORT as BCFTOOLS_SORT2 } from '../../../modules/nf-core/bcftools/sort'
include { TABIX_TABIX as TABIX_TABIX_1    } from '../../../modules/nf-core/tabix/tabix'
include { TABIX_TABIX as TABIX_TABIX_2    } from '../../../modules/nf-core/tabix/tabix'


workflow SV_VCF_CONVERSIONS {
    take:
    input_ch    // channel: [val(meta), vcf]
    fai         // reference channel [val(meta), ref.fa.fai]

    main:

    if (params.sv_standardization.contains("variantextractor")){
        // uses VariantExtractor to homogenize variants
        VARIANTEXTRACTOR(
            input_ch
        )

        // sort vcf
        BCFTOOLS_SORT1(
            VARIANTEXTRACTOR.out.vcf
        )
        input_ch = BCFTOOLS_SORT1.out.vcf
    }

    if (params.sv_standardization.contains("svtk")){

        out_vcf_ch = channel.empty()

        supported_callers2 = ["delly", "melt", "manta", "wham", "dragen", "lumpy", "scrable", "smoove"]
        input_ch
            .branch{ meta, _vcf->
                def caller = meta.caller
                def supported = supported_callers2.contains(caller)
                if(!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported in svtk. Skipping standardization...")
                }
                tool:  supported
                other: !supported
            }
            .set{input}

        TABIX_TABIX_1(
            input.tool
        )

        SVTK_STANDARDIZE(
            input.tool.join(TABIX_TABIX_1.out.index),
            fai
        )

        BCFTOOLS_SORT2(
            SVTK_STANDARDIZE.out.vcf
        )

        out_vcf_ch.mix(
                BCFTOOLS_SORT2.out.vcf,
                input.other
            ).set{input_ch}

    }

    if (params.sv_standardization.contains("svdecompose")){
        RTGTOOLS_SVDECOMPOSE(
            input_ch.map{ meta, vcf -> tuple(meta, vcf, [])}
        )
        input_ch = RTGTOOLS_SVDECOMPOSE.out.vcf
    }

    input_ch
        .branch { files ->
            compressed:   files[1].getName().endsWith('.gz')
            uncompressed: true
        }
        .set { ch_inputs }

    TABIX_BGZIPTABIX (
        ch_inputs.uncompressed
        )

    TABIX_TABIX_2 (
        ch_inputs.compressed
        )

    compressed_ch = ch_inputs.compressed.join(TABIX_TABIX_2.out.index)
    vcf_ch = TABIX_BGZIPTABIX.out.gz_index.mix(compressed_ch)

    // RUN SVYNC tool to reformat SV callers
    if(params.sv_standardization.contains("svync")){
        out_vcf_ch = channel.empty()
        supported_callers = ["delly", "dragen", "gridss", "manta", "smoove"]

        vcf_ch
            .branch{ meta, vcf, tbi ->
                def caller = meta.caller
                def supported = supported_callers.contains(caller)
                if(!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported in svync. Skipping standardization...")
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
        out_vcf_ch.mix(
                SVYNC.out.vcf,
                input.other
            )
            .map{ file ->
                def meta = file[0]
                def vcf = file[1]
                [ meta, vcf ]
            }
            .set { vcf_ch }
    }



    emit:
    vcf_ch   // channel: [val(meta), vcf]
}
