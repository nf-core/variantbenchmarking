//
// CONCORDANCE_ANALYSIS: SUBWORKFLOW FOR CONCORDANCE ANALYSIS BETWEEN BASE AND TEST VCFS
//

include { GATK4_CONCORDANCE                      } from '../../../modules/nf-core/gatk4/concordance'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_TP_BASE } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FN      } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_FP      } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_TP_COMP } from '../../../modules/nf-core/bcftools/view'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_TP_BASE } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_FN      } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_FP      } from '../../../modules/nf-core/bcftools/reheader'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_TP_COMP } from '../../../modules/nf-core/bcftools/reheader'

workflow CONCORDANCE_ANALYSIS {
    take:
    input_ch   // channel: [val(meta), vcf.gz, index]
    bed_ch     // channel: [path(bed)] 
    fasta_ch   // reference channel [val(meta), ref.fa]  
    fai_ch     // reference channel [val(meta), ref.fa.fai]
    dictionary // reference channel [val(meta), ref.dict]

    main:

    tagged_variants = channel.empty()

    ch_pairs = input_ch
        .map { meta, vcf1, tbi1 ->
            // Simplify to just keep the id and files
            [meta.id, vcf1, tbi1]
        }
        .toList()
        .flatMap { items ->
            def result = []

            // Generate pairwise combinations
            for (int i = 0; i < items.size(); i++) {
                for (int j = i + 1; j < items.size(); j++) {
                    def left = items[i]   // [id1, vcf1, tbi1]
                    def right = items[j]  // [id2, vcf2, tbi2]

                    // Create new metadata with combined IDs
                    def combinedMeta = [id: "${left[0]}-${right[0]}"]

                    result << [
                        combinedMeta,     // [id: "test7-test6"]
                        left[1], left[2], // vcf1, tbi1 from first sample
                        right[1], right[2] // vcf2, tbi2 from second sample
                    ]
                }
            }
            return result
        }

    ch_bed_input = bed_ch
        .map { file -> tuple(["id": "intervals"], file) }
        .ifEmpty([[:], []])

    // GATK4 concordance does not support structural variants now - GATK4 SVCONCORDANCE is in beta
    GATK4_CONCORDANCE(
        ch_pairs,
        ch_bed_input,
        fasta_ch,
        fai_ch,
        dictionary
    )

    // tag meta and collect summary reports
    GATK4_CONCORDANCE.out.summary
        .map { _meta, file -> tuple([vartype: params.variant_type] + [benchmark_tool: "concordance"], file) }
        .groupTuple()
        .set{ summary_reports }

    // Split FN variants from TPFN variants
    BCFTOOLS_VIEW_FN(
        GATK4_CONCORDANCE.out.tpfn.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []
    )

    // Reheader FN variants
    BCFTOOLS_REHEADER_FN(

        BCFTOOLS_VIEW_FN.out.vcf.map{ meta, file ->
            [ meta, file, [], [] ]
        },
        fai_ch
    )

    // Tag FN variants
    BCFTOOLS_REHEADER_FN.out.vcf
        .join(BCFTOOLS_REHEADER_FN.out.index)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FN"] + [id: "concordance"], file, index) }
        .set { vcf_fn }

    // Split TP base variants from TP-FN variants
    BCFTOOLS_VIEW_TP_BASE(
        GATK4_CONCORDANCE.out.tpfn.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []
    )

    // Reheader TP variants
    BCFTOOLS_REHEADER_TP_BASE(

        BCFTOOLS_VIEW_TP_BASE.out.vcf.map{ meta, file ->
            [ meta, file, [], [] ]
        },
        fai_ch
    )

    // Tag TP variants
    BCFTOOLS_REHEADER_TP_BASE.out.vcf
        .join(BCFTOOLS_REHEADER_TP_BASE.out.index)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_base"] + [id: "concordance"], file, index) }
        .set { vcf_tp_base }

    // Split TP comp variants from TPFP variants
    BCFTOOLS_VIEW_TP_COMP(
        GATK4_CONCORDANCE.out.tpfp.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []
    )

    // Reheader TP comp variants
    BCFTOOLS_REHEADER_TP_COMP(

        BCFTOOLS_VIEW_TP_COMP.out.vcf.map{ meta, file ->
            [ meta, file, [], [] ]
        },
        fai_ch
    )

    // Tag TP comp variants
    BCFTOOLS_REHEADER_TP_COMP.out.vcf
        .join(BCFTOOLS_REHEADER_TP_COMP.out.index)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "TP_comp"] + [id: "concordance"], file, index) }
        .set { vcf_tp_comp }

    // Split FP variants from TPFP variants
    BCFTOOLS_VIEW_FP(
        GATK4_CONCORDANCE.out.tpfp.map{ meta, vcf -> tuple(meta, vcf, []) },
        [],
        [],
        []
    )

    // Reheader FP variants
    BCFTOOLS_REHEADER_FP(
        BCFTOOLS_VIEW_FP.out.vcf.map{ meta, file ->
            [ meta, file, [], [] ]
        },
        fai_ch
    )

    // Tag FP variants
    BCFTOOLS_REHEADER_FP.out.vcf
        .join(BCFTOOLS_REHEADER_FP.out.index)
        .map { _meta, file, index -> tuple([vartype: params.variant_type] + [tag: "FP"] + [id: "concordance"], file, index) }
        .set { vcf_fp }

    tagged_variants = tagged_variants.mix(
        vcf_fn,
        vcf_fp,
        vcf_tp_base,
        vcf_tp_comp
    )

    emit:
    summary_reports // channel: [val(meta), reports]
    tagged_variants // channel: [val(meta), vcfs]
}
