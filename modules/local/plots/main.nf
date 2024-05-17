process PLOTS {
    tag "$benchmark_tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-4efdedce78cbbcee99652df1954093caf4730f6a:6da101c870eae13e6ddbf6940d3797d9b86ae763':
        'biocontainers/mulled-v2-4efdedce78cbbcee99652df1954093caf4730f6a:6da101c870eae13e6ddbf6940d3797d9b86ae763' }"

    input:
    tuple val(benchmark_tool), path(summary)

    output:
    tuple val(benchmark_tool),path("*.png"), emit: plots
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    plots.R $summary $benchmark_tool

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
