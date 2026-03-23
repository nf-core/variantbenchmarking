
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge TP/FP/FN results from different tools.
//

include { GAWK as REFORMAT_HEADER          } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIP as TABIX_BGZIP_UNZIP } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_BGZIPTABIX                 } from '../../../modules/nf-core/tabix/bgziptabix'
include { BCFTOOLS_MERGE                   } from '../../../modules/nf-core/bcftools/merge'
include { SURVIVOR_MERGE                   } from '../../../modules/nf-core/survivor/merge'
include { VCF_TO_CSV                       } from '../../../modules/local/custom/vcf_to_csv'
include { SOMPY_FEATURES_MERGE             } from '../../../modules/local/sompy_features/merge'
include { PLOTS_UPSET                      } from '../../../modules/local/plots/upset'


workflow COMPARE_BENCHMARK_RESULTS {
    take:
    evaluations     // channel: [val(meta), vcf.gz, index]
    evaluations_csv // channel: [val(meta), csv]
    fasta           // reference channel [val(meta), ref.fa]
    fai             // reference channel [val(meta), ref.fa.fai]

    main:
    merged_vcfs = channel.empty()
    ch_plots    = channel.empty()

    if (params.variant_type == "small" | params.variant_type == "snv" | params.variant_type == "indel"){

        // Small Variants
        REFORMAT_HEADER(
            evaluations.map { meta, vcf, _tbi -> [meta, vcf] },
            [],
            false
        )

        TABIX_BGZIPTABIX(
            REFORMAT_HEADER.out.output
        )

        // merge small variants
        BCFTOOLS_MERGE(
            TABIX_BGZIPTABIX.out.gz_index.groupTuple(),
            fasta,
            fai,
            [[],[]]
        )
        merged_vcfs = merged_vcfs.mix(BCFTOOLS_MERGE.out.vcf)
    }
    else{
        // SV part
        // unzip vcfs
        TABIX_BGZIP_UNZIP(
            evaluations.map { item -> tuple(item[0], item[1]) }
        )

        TABIX_BGZIP_UNZIP.out.output
            .groupTuple()
            .set{vcf_ch}

        // Merge Benchmark SVs from different tools
        SURVIVOR_MERGE(
            vcf_ch,
            1000,
            1,
            1,
            0,
            0,
            30
        )
        merged_vcfs = merged_vcfs.mix(SURVIVOR_MERGE.out.vcf)

    }

    // convert vcf files to csv
    VCF_TO_CSV(
        merged_vcfs
    )

    SOMPY_FEATURES_MERGE(
        evaluations_csv.groupTuple()
    )

    if (!params.skip_plots.contains("upset")){
        VCF_TO_CSV.out.output.mix(SOMPY_FEATURES_MERGE.out.output).map{
            meta, csv ->
                def newMeta = meta.clone()
                newMeta.remove('tag')
            tuple(newMeta,csv)
        }.set{upset_input}

        PLOTS_UPSET(
            upset_input.groupTuple()
        )
        ch_plots = ch_plots.mix(PLOTS_UPSET.out.plot)
    }

    emit:
    merged_vcfs  // channel: [val(meta), vcf]
    ch_plots     // channel: [.png]

}
