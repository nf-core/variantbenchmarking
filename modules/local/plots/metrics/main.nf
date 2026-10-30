process PLOTS_METRICS {
    tag "$meta.benchmark_tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/577648e514e596e9c07da1d91ccc7c60ad63fcbfe5b658aee73d3adca9cd337b/data' :
        'community.wave.seqera.io/library/python_pip_pandas_plotly_upsetplot:1131451904bdd81f' }"

    input:
    tuple val(meta), path(summary)

    output:
    path("*.html") , emit: plots
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"

    """
    plot_metrics.py \\
        $summary \\
        ${meta.benchmark_tool} \\
        --output variant_metrics_${prefix}.html \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"
    """
    touch variant_metrics_${prefix}.html
    """
}
