//
// PREPARE_VCFS: SUBWORKFLOW TO PREPARE INPUT VCFS
//

params.options = [:]

include { BGZIP_TABIX      } from '../../modules/local/bgzip_tabix.nf'       addParams( options: params.options )
include { BCFTOOLS_VIEW    } from '../../modules/local/bcftools_view'      addParams( options: params.options )
include { BCFTOOLS_NORM as BCFTOOLS_NORM_1 } from '../../modules/nf-core/bcftools/norm'      addParams( options: params.options )
include { BCFTOOLS_NORM as BCFTOOLS_NORM_2 } from '../../modules/nf-core/bcftools/norm'      addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_1   } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_2   } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { TABIX_TABIX   as TABIX_TABIX_3   } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { BGZIP_TABIX as BGZIP_TABIX_1     } from '../../modules/local/bgzip_tabix'      addParams( options: params.options )
include { BGZIP_TABIX as BGZIP_TABIX_2     } from '../../modules/local/bgzip_tabix'      addParams( options: params.options )

workflow PREPARE_VCFS_TRUTH {
    take:
    truth_ch    // channel: [val(meta), vcf]
    ref         // reference channel [ref.fa, ref.fa.fai]
    main_chroms // channel: path(chrom sizes)

    main:

    versions=Channel.empty()

// Check tool spesific conversions

    // https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
    // GRIDSS simple event annotation

    emit:
    vcf_ch
    versions
}
