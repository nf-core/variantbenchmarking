process PLOTS {
    tag "$benchmark_tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-400bee8f069b1e20921115d7a92e5e6dd83bc49c:055769a32b18e19ed8bee8ab152eedac16879359-0':
        'biocontainers/mulled-v2-400bee8f069b1e20921115d7a92e5e6dd83bc49c:055769a32b18e19ed8bee8ab152eedac16879359-0' }"

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
