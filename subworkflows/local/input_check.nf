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
    meta.id           = row.caller
    meta.vartype      = row.vartype

    // add path(s) of the fastq file(s) to the meta map
    def vcf_meta = []
        if (!file(row.test_vcf).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Test file does not exist!\n${row.test_vcf}"
        }
        if (meta.vartype == "sv"){
            if (meta.id == "delly"){
                vcf_meta = [  meta, file(row.test_vcf), file("${projectDir}/assets/svync/delly.yaml")]
            }
            else if (meta.id == "gridss"){
                vcf_meta = [  meta,  file(row.test_vcf), file("${projectDir}/assets/svync/gridss.yaml")]
            }
            else if (meta.id == "manta"){
                if (file("${projectDir}/assets/svync/manta.yaml").exists()){
                    vcf_meta = [  meta,  file(row.test_vcf), file("${projectDir}/assets/svync/manta.yaml")]
                }
            }
            else if (meta.id == "smoove"){
                vcf_meta = [  meta, file(row.test_vcf), file("${projectDir}/assets/svync/smoove.yaml")]
            }
            else{
                vcf_meta = [  meta, file(row.test_vcf), file("${projectDir}/assets/svync/default.yaml")]
            }
        }else{
            vcf_meta = [  meta, file(row.test_vcf)]
        }

    return vcf_meta
}
