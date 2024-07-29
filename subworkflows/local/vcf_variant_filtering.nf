//
// VCF_VARIANT_FILTERING: Filter SV using survivor and bcftools
//

include { SURVIVOR_FILTER     } from '../../modules/nf-core/survivor/filter'
include { TABIX_BGZIP         } from '../../modules/nf-core/tabix/bgzip'
include { TABIX_BGZIPTABIX    } from '../../modules/nf-core/tabix/bgziptabix'
include { BCFTOOLS_FILTER     } from '../../modules/nf-core/bcftools/filter'


workflow VCF_VARIANT_FILTERING {
    take:
    vcf_ch    // channel: [val(meta),vcf.gz, tbi]

    main:

    versions = Channel.empty()

    //
    // TABIX_BGZIP
    //
    // unzip vcf file, required for survivor filter
    TABIX_BGZIP(
        vcf_ch.map { meta, vcf, tbi -> [ meta, vcf ]}
    )
    versions = versions.mix(TABIX_BGZIP.out.versions.first())
    vcf_ch = TABIX_BGZIP.out.output

    if(params.exclude_expression  != null & params.include_expression  != null){
        //
        // BCFTOOLS_FILTER
        //
        // filter vcf files using bcftools expressions
        BCFTOOLS_FILTER(
            vcf_ch
        )
        versions = versions.mix(BCFTOOLS_FILTER.out.versions.first())
        vcf_ch = BCFTOOLS_FILTER.out.vcf
    }

    out_vcf_ch = Channel.empty()

    if(params.min_sv_size > 0 | params.max_sv_size != -1 | params.min_allele_freq != -1 | params.min_num_reads != -1){
        vcf_ch.branch{
                def meta = it[0]
                sv:  meta.vartype == "sv"
                small:  meta.vartype == "small"
                cnv:  meta.vartype == "cnv"
                snv: meta.vartype == "snv"
                indel: meta.vartype == "indel"
                other: false
            }
            .set{vcf}

        //
        // MODULE: SURVIVOR_FILTER
        //
        // filters out smaller SVs than min_sv_size, only applicable to SV files
        SURVIVOR_FILTER(
            vcf.sv.map{ meta, vcf -> 
                [ meta, vcf, [] ]
            },
            params.min_sv_size,
            params.max_sv_size,
            params.min_allele_freq,
            params.min_num_reads
        )
        versions = versions.mix(SURVIVOR_FILTER.out.versions.first())

        out_vcf_ch = out_vcf_ch.mix(
            SURVIVOR_FILTER.out.vcf,
            vcf.small,
            vcf.cnv,
            vcf.snv,
            vcf.indel
        )
        vcf_ch = out_vcf_ch
    }
    //
    // TABIX_BGZIPTABIX
    //
    // zip and index vcf files
    TABIX_BGZIPTABIX(
        vcf_ch
    )
    versions = versions.mix(TABIX_BGZIPTABIX.out.versions.first())
    vcf_ch = TABIX_BGZIPTABIX.out.gz_tbi

    emit:
    vcf_ch
    versions
}
