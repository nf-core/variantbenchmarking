//
// LIFTOVER_VCFS: SUBWORKFLOW TO LIFTOVER VCFS HG37 TO HG38 OR HG38 TO HG37
//

include { PICARD_LIFTOVERVCF   } from '../../../modules/nf-core/picard/liftovervcf'
include { REFORMAT_HEADER      } from '../../../modules/local/custom/reformat_header'
include { BCFTOOLS_ANNOTATE    } from '../../../modules/nf-core/bcftools/annotate'
include { UCSC_LIFTOVER        } from '../../../modules/nf-core/ucsc/liftover'
include { GNU_SORT             } from '../../../modules/nf-core/gnu/sort'
include { BEDTOOLS_MERGE       } from '../../../modules/nf-core/bedtools/merge'


workflow LIFTOVER_VCFS {
    take:
    ch_vcf          // channel: [val(meta), vcf]
    ch_bed          // channel: [bed]
    fasta           // reference channel [val(meta), ref.fa]
    chain           // chain channel [val(meta), chain.gz]
    rename_chr      // reference channel [val(meta), chrlist.txt]
    dictionary      // reference channel [val(meta), genome.dict]

    main:

    versions = Channel.empty()

    // Use picard liftovervcf tool to convert vcfs
    PICARD_LIFTOVERVCF(
        ch_vcf,
        dictionary,
        fasta,
        chain
    )
    versions = versions.mix(PICARD_LIFTOVERVCF.out.versions)
    vcf_ch   = PICARD_LIFTOVERVCF.out.vcf_lifted

    // reformat header, convert PS TYPE integer to string after liftover
    REFORMAT_HEADER(
        vcf_ch.map{meta, vcf -> tuple(meta, vcf, [])}
    )
    versions = versions.mix(REFORMAT_HEADER.out.versions)

    // rename chr after liftover
    BCFTOOLS_ANNOTATE(
        REFORMAT_HEADER.out.gz_tbi.map{meta, vcf, tbi -> tuple(meta, vcf, tbi, [], [])},
        [],
        [],
        rename_chr.map{_meta, file -> file}
    )
    vcf_ch = BCFTOOLS_ANNOTATE.out.vcf

    // liftover high confidence bed file if given
    UCSC_LIFTOVER(
        ch_bed.map{file -> tuple([id: params.truth_id], file)},
        chain.map{_meta, file -> file}
    )
    versions = versions.mix(UCSC_LIFTOVER.out.versions)

    // sort bed file
    GNU_SORT(
        UCSC_LIFTOVER.out.lifted
    )

    // merge the intersected regions
    BEDTOOLS_MERGE(
        GNU_SORT.out.bed
    )
    versions = versions.mix(BEDTOOLS_MERGE.out.versions)
    bed_ch = BEDTOOLS_MERGE.out.bed

    emit:
    vcf_ch      // channel: [val(meta), vcf.gz]
    bed_ch      // channel: [val(meta), bed]
    versions    // channel: [versions.yml]
}
