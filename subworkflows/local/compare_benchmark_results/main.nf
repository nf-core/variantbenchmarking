
//
// COMPARE_BENCHMARK_RESULTS: SUBWORKFLOW to merge TP/FP/FN results from different tools.
//

include { GAWK as REFORMAT_HEADER                  } from '../../../modules/nf-core/gawk'
include { TABIX_BGZIP as TABIX_BGZIP_UNZIP         } from '../../../modules/nf-core/tabix/bgzip'
include { TABIX_BGZIPTABIX                         } from '../../../modules/nf-core/tabix/bgziptabix'
include { BCFTOOLS_MERGE                           } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_INDEX                           } from '../../../modules/nf-core/bcftools/index'
include { SURVIVOR_MERGE                           } from '../../../modules/nf-core/survivor/merge'
include { GATK4_VARIANTSTOTABLE as VARIANTSTOTABLE } from '../../../modules/nf-core/gatk4/variantstotable'
include { MERGE_SOMPY_FEATURES                     } from '../../../modules/local/custom/merge_sompy_features'
include { PLOT_UPSET                               } from '../../../modules/local/custom/plot_upset'


workflow COMPARE_BENCHMARK_RESULTS {
    take:
    evaluations     // channel: [val(meta), vcf.gz, index]
    evaluations_csv // channel: [val(meta), csv]
    fasta           // reference channel [val(meta), ref.fa]
    fai             // reference channel [val(meta), ref.fa.fai]
    dictionary      // reference channel [val(meta), genome.dict]

    main:
    versions    = channel.empty()
    merged_vcfs = channel.empty()
    merged_tbis = channel.empty()
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
        versions = versions.mix(BCFTOOLS_MERGE.out.versions_bcftools.first())

        merged_vcfs = merged_vcfs.mix(BCFTOOLS_MERGE.out.vcf)
        merged_tbis = merged_tbis.mix(BCFTOOLS_MERGE.out.index)
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

        // index merged vcf file
        BCFTOOLS_INDEX(
        merged_vcfs
        )
        versions = versions.mix(BCFTOOLS_INDEX.out.versions_bcftools.first())

        merged_tbis = merged_tbis.mix(BCFTOOLS_INDEX.out.tbi)
    }

    variantstotable_input_ch = merged_vcfs
                                .join(merged_tbis)
                                .map{ meta, vcf, tbi -> [meta, vcf, tbi, [], [], []] }

    // convert vcf files to tsv
    VARIANTSTOTABLE(
    variantstotable_input_ch,
    fasta,
    fai,
    dictionary
    )
    versions = versions.mix(VARIANTSTOTABLE.out.versions_gatk4.first())

    MERGE_SOMPY_FEATURES(
        evaluations_csv.groupTuple()
    )
    versions = versions.mix(MERGE_SOMPY_FEATURES.out.versions.first())

    if (!params.skip_plots.contains("upset")){
        VARIANTSTOTABLE.out.table.mix(MERGE_SOMPY_FEATURES.out.output).map{
            meta, tsv ->
                def newMeta = meta.clone()
                newMeta.remove('tag')
            tuple(newMeta,tsv)
        }.set{upset_input}

        PLOT_UPSET(
            upset_input.groupTuple()
        )
        versions = versions.mix(PLOT_UPSET.out.versions)
        ch_plots = ch_plots.mix(PLOT_UPSET.out.plot)
    }

    emit:
    merged_vcfs  // channel: [val(meta), vcf]
    ch_plots     // channel: [.png]
    versions     // channel: [versions.yml]

}
