process TRUVARI_REFINE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/truvari:4.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/truvari:4.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), val(meta2), path(bench)
    each path(bed)
    tuple path(fasta), path(fai)

    output:
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    truvari refine \\
        --use-original-vcfs \\
        --reference $fasta \\ 
        --regions $bed \\ 
        $bench
        . 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        truvari: \$(echo \$(truvari version 2>&1) | sed 's/^Truvari v//' ))
    END_VERSIONS
    """
}
