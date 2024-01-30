process TRUVARI_BENCH {
    tag "$meta.id $meta2.caller"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/truvari:4.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/truvari:4.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta),val(meta2), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi)
    each path(bed)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta),val(meta2), path(bench)  , emit: bench
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--includebed $bed" : ""

    """
    truvari bench \\
        --base ${truth_vcf} \\
        --comp ${vcf} \\
        --reference ${fasta} \\
        --output ${prefix} \\
        ${regions} \\
        ${args}

    mkdir bench
    mv ${prefix}/fn.vcf.gz          bench/
    mv ${prefix}/fn.vcf.gz.tbi      bench/
    mv ${prefix}/fp.vcf.gz          bench/
    mv ${prefix}/fp.vcf.gz.tbi      bench/
    mv ${prefix}/tp-base.vcf.gz     bench/
    mv ${prefix}/tp-base.vcf.gz.tbi bench/
    mv ${prefix}/tp-comp.vcf.gz     bench/
    mv ${prefix}/tp-comp.vcf.gz.tbi bench/
    mv ${prefix}/summary.json       bench/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' ))
    END_VERSIONS
    """
}
