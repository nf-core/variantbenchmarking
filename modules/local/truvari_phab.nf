process TRUVARI_PHAB {
    tag "$meta.id $meta2.caller"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/truvari:v4.3.0':
        'kubran/truvari:v4.3.0' }"

    input:
    tuple val(meta),  path(vcf), path(tbi), path(truth_vcf), path(truth_tbi), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta),path("*.vcf.gz")         , emit: harmon
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--region $bed" : ""

    """
    truvari phab \\
        --base ${truth_vcf} \\
        --comp ${vcf} \\
        --bSample $meta.id \\
        --cSample $meta.id \\
        --reference ${fasta} \\
        --output ${prefix}.vcf.gz \\
        ${regions} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' ))
    END_VERSIONS
    """
}
