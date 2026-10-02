//
// SV_VCF_CONVERSIONS: SUBWORKFLOW to apply tool spesific conversions
//

include { SVYNC                           } from '../../../modules/nf-core/svync'
include { VARIANTEXTRACTOR                } from '../../../modules/nf-core/variantextractor'
include { SVTK_STANDARDIZE                } from '../../../modules/nf-core/svtk/standardize'
include { RTGTOOLS_SVDECOMPOSE            } from '../../../modules/nf-core/rtgtools/svdecompose'
include { BCFTOOLS_SORT as BCFTOOLS_SORT1 } from '../../../modules/nf-core/bcftools/sort'
include { BCFTOOLS_SORT as BCFTOOLS_SORT2 } from '../../../modules/nf-core/bcftools/sort'


workflow SV_VCF_CONVERSIONS {
    take:
    input_ch // channel: [val(meta), vcf]
    fai // reference channel [val(meta), ref.fa.fai]

    main:

    vcf_ch = input_ch

    if (params.sv_standardization.contains("variantextractor")) {
        // uses VariantExtractor to homogenize variants
        VARIANTEXTRACTOR(
            vcf_ch
        )

        // sort vcf
        BCFTOOLS_SORT1(
            VARIANTEXTRACTOR.out.vcf
        )
        vcf_ch = BCFTOOLS_SORT1.out.vcf
    }

    if (params.sv_standardization.contains("svtk")) {

        svtk_callers = ["delly", "melt", "manta", "wham", "dragen", "lumpy", "scrable", "smoove"]
        vcf_ch
            .branch { meta, _vcf ->
                def caller = meta.caller
                def supported = svtk_callers.contains(caller)
                if (!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported in svtk. Skipping standardization...")
                }
                supported: supported
                unsupported: !supported
            }
            .set { svtk_input }

        SVTK_STANDARDIZE(
            svtk_input.supported,
            fai,
        )

        BCFTOOLS_SORT2(
            SVTK_STANDARDIZE.out.vcf
        )

        vcf_ch = BCFTOOLS_SORT2.out.vcf.mix(svtk_input.unsupported)
    }

    if (params.sv_standardization.contains("svdecompose")) {
        RTGTOOLS_SVDECOMPOSE(
            vcf_ch.map { meta, vcf -> tuple(meta, vcf, []) }
        )
        vcf_ch = RTGTOOLS_SVDECOMPOSE.out.vcf
    }

    // RUN SVYNC tool to reformat SV callers
    if (params.sv_standardization.contains("svync")) {
        svync_callers = ["delly", "dragen", "gridss", "manta", "smoove"]

        vcf_ch
            .branch { meta, _vcf ->
                def caller = meta.caller
                def supported = svync_callers.contains(caller)
                if (!supported) {
                    log.warn("Standardization for SV caller '${caller}' is not supported in svync. Skipping standardization...")
                }
                supported: supported
                unsupported: !supported
            }
            .set { svync_input }

        SVYNC(
            svync_input.supported.map { meta, vcf ->
                [meta, vcf, [], file("${projectDir}/assets/svync/${meta.caller}.yaml", checkIfExists: true)]
            }
        )

        vcf_ch = SVYNC.out.vcf.mix(svync_input.unsupported)
    }

    emit:
    vcf_ch // channel: [val(meta), vcf]
}
