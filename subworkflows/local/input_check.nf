//
// Check input samplesheet and get read channels
//

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK (samplesheet)
        .csv
        .splitCsv ( header:true, sep:',' )
        .map{ create_vcf_channel(it) }
        .set {ch_sample}

    emit:
    ch_sample // channel: [ val(meta), test_vcf]
    versions = SAMPLESHEET_CHECK.out.versions
}

def create_vcf_channel(LinkedHashMap row) {
// create meta map
    def meta = [:]
    meta.id           = params.sample

    def meta2 = [:]
    meta2.caller      = row.caller
    // add path(s) of the fastq file(s) to the meta map
    def vcf_meta = []
        if (!file(row.test_vcf).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Test file does not exist!\n${row.test_vcf}"
        }
        vcf_meta = [  meta, meta2, file(row.test_vcf)]
    return vcf_meta
}
