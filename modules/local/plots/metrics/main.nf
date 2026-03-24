process PLOTS_METRICS {
    tag "$meta.benchmark_tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs:1.8.4--r43hdfd78af_0':
        'biocontainers/r-shinyngs:1.8.4--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(summary)

    output:
    path("*.png")          , emit: plots
    tuple val("${task.process}"), val('r_base'), eval("R --version 2>&1 | head -n 1 | sed 's/^R version //; s/ .*//'"), topic: versions, emit: versions_rbase

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"

    """
    plots.R $summary $meta.benchmark_tool $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"
    """
    touch metric_by_tool_${prefix}.png
    touch variants_by_tool_${prefix}.png
    """
}
