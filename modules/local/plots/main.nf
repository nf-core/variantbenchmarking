process PLOTS {
    tag "$meta.benchmark_tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs:1.8.4--r43hdfd78af_0':
        'biocontainers/r-shinyngs:1.8.4--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(summary)

    output:
    tuple val(meta),path("*.png"), emit: plots
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"

    """
    plots.R $summary $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"
    """
    touch metric_by_tool_${prefix}.png
    touch variants_by_tool_${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
