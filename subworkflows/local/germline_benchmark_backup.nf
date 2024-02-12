//
// GERMLINE: SUBWORKFLOW FOR GERMLINE VARIANTS
//

params.options = [:]

include { TRUVARI_BENCH          } from '../../modules/nf-core/truvari/bench'          addParams( options: params.options )
include { SVANALYZER_SVBENCHMARK } from '../../modules/nf-core/svanalyzer/svbenchmark' addParams( options: params.options )
include { WITTYER                } from '../../modules/local/wittyer'                  addParams( options: params.options )
include { VCFDIST                } from '../../modules/local/vcfdist'                  addParams( options: params.options )
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_QUERY          } from '../../modules/local/bcftools_view'       addParams( options: params.options )
include { BCFTOOLS_VIEW as BCFTOOLS_VIEW_TRUTH          } from '../../modules/local/bcftools_view'       addParams( options: params.options )
include { ADDHEAD as ADDHEAD_TRUTH  } from '../../modules/local/addhead' addParams( options: params.options )
include { ADDHEAD as ADDHEAD_QUERY  } from '../../modules/local/addhead' addParams( options: params.options )
include { BCFTOOLS_ISEC as BCFTOOLS_ISEC_TRUTH          } from '../../modules/nf-core/bcftools/isec'     addParams( options: params.options )
include { BCFTOOLS_ISEC as BCFTOOLS_ISEC_QUERY          } from '../../modules/nf-core/bcftools/isec'     addParams( options: params.options )
include { TABIX_TABIX as TABIX_TABIX_1  } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )
include { TABIX_TABIX as TABIX_TABIX_2  } from '../../modules/nf-core/tabix/tabix'        addParams( options: params.options )

workflow GERMLINE_BENCHMARK {
    take:
    input_ch  // channel: [val(meta),val(meta2), test_vcf, test_index , truth_vcf, truth_index]
    bed       // channel: bed
    ref       // reference channel [ref.fa, ref.fa.fai]
    truth_vcf // channel: [val(meta),val(meta2),truth_vcf, truth_index]

    main:

    versions=Channel.empty()

    // SV benchmarking

    if (params.method.contains('truvari')){
        //
        // MODULE: TRUVARI_BENCH
        //
        TRUVARI_BENCH(
            input_ch,
            bed,
            ref
        )
        versions = versions.mix(TRUVARI_BENCH.out.versions)
    }

    if (params.method.contains('svanalyzer')){
        //
        // MODULE: SVANALYZER_SVBENCHMARK
        //
        // slower than truvari
        SVANALYZER_SVBENCHMARK(
            input_ch,
            ref,
            bed        
        )
        versions = versions.mix(SVANALYZER_SVBENCHMARK.out.versions)

        // get the original headers from the vcfs
        //
        // MODULE: BCFTOOLS_VIEW
        //
        BCFTOOLS_VIEW_QUERY(
            input_ch.map{it -> tuple(it[0],it[1], it[2], it[3])}
        )
        query_header = BCFTOOLS_VIEW_QUERY.out.header

        BCFTOOLS_VIEW_TRUTH(
            truth_vcf
        )
        truth_header = BCFTOOLS_VIEW_TRUTH.out.header
        versions = versions.mix(BCFTOOLS_VIEW_QUERY.out.versions)
        //
        // MODULE: BCFTOOLS_REHEADER
        //
        SVANALYZER_SVBENCHMARK.out.fns.combine(truth_header, by:0)
                                        .map{it -> tuple(it[0],it[1], it[2], it[4])}
                                        .set{fns_header}
        fns_header.view()
        ADDHEAD_TRUTH(
            fns_header
        )

        TABIX_TABIX_1(
            ADDHEAD_TRUTH.out.vcf
        )

        ADDHEAD_TRUTH.out.vcf.join(TABIX_TABIX_1.out.tbi, by:1)
                        .map{it -> tuple( it[1], it[0], it[2], it[4])}
                        .set{fns_vcf}

        //////
        SVANALYZER_SVBENCHMARK.out.fps.combine(query_header, by:0)
                .map{it -> tuple(it[0],it[1], it[2], it[4])}
                .set{fps_header}

        ADDHEAD_QUERY(
            fps_header
        )
        versions = versions.mix(ADDHEAD_QUERY.out.versions)

        TABIX_TABIX_2(
            ADDHEAD_QUERY.out.vcf
        )

        ADDHEAD_QUERY.out.vcf.join(TABIX_TABIX_2.out.tbi, by:1)
                        .map{it -> tuple( it[1], it[0], it[2], it[4])}
                        .set{fps_ch}
        //
        // MODULE: BCFTOOLS_ISEC
        // 

        // Find TP_comp (query)
        input_ch.map{it -> tuple(it[0],it[1], it[2], it[3])}
                .combine(fps_ch, by:0)
                .map{it -> tuple(it[0],it[1], it[2], it[5], it[3], it[6])}
                .set{query_ch}

        BCFTOOLS_ISEC_QUERY(
            query_ch
        )

        // Find TP_comp (query)
        truth_vcf.combine(fns_vcf, by:0)
                .map{it -> tuple(it[0],it[1], it[2], it[5], it[3], it[6])}
                .set{truth_ch}

        // Find TP_base (truth)
        BCFTOOLS_ISEC_TRUTH(
            truth_ch
        )
        versions = versions.mix(BCFTOOLS_ISEC_TRUTH.out.versions)

    } 

    if (params.method.contains('wittyer')){
        //
        // MODULE: WITTYER
        //
        WITTYER(
            input_ch,
            bed,
            []
        )
        versions = versions.mix(WITTYER.out.versions)
    }


    if (params.method.contains('vcfdist')){
        //
        // MODULE: VCFDIST
        //
        VCFDIST(
            input_ch,
            ref,
            bed
        )
        versions = versions.mix(VCFDIST.out.versions)
    }


    emit:
    versions
}
