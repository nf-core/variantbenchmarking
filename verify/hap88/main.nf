nextflow.enable.dsl = 2

process HAP88_VCF_TO_CSV {
    tag "${prefix}"
    container 'community.wave.seqera.io/library/pip_pandas:40d2e76c16c136f0'
    cpus 1
    memory '6.GB'
    errorStrategy 'finish'
    maxRetries 0

    input:
    tuple val(prefix), path(vcf)
    path(converter)

    output:
    path("${prefix}.csv")

    script:
    """
    python3 - << 'HASHCHECK'
import hashlib, pathlib, sys
p = pathlib.Path("${converter}")
data = p.read_bytes()
header = 'blob %d' % len(data)
blob = header.encode() + bytes([0]) + data
digest = hashlib.sha1(blob).hexdigest()
print('vcf_to_csv.py blob', digest, 'bytes', len(data))
if digest != '28619cb7e23d7add57efd70f5eafb1fcf08949c5':
    sys.exit(2)
HASHCHECK
    if command -v python3 >/dev/null 2>&1; then PYBIN=python3; else PYBIN=python; fi
    "\$PYBIN" "${converter}" "${vcf}" "${prefix}.csv"
    ls -l "${prefix}.csv"
    """
}

workflow {
    converter = file("${projectDir}/bin/vcf_to_csv.py", checkIfExists: true)
    Channel
        .of(
            ['happy.TP_base', params.base_vcf],
            ['happy.TP_comp', params.comp_vcf]
        )
        .map { prefix, uri -> tuple(prefix, file(uri, checkIfExists: true)) }
        | set { ch }
    HAP88_VCF_TO_CSV(ch, converter)
}
