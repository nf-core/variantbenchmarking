process PLOTS_SVLEN_DIST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/57/577648e514e596e9c07da1d91ccc7c60ad63fcbfe5b658aee73d3adca9cd337b/data' :
        'community.wave.seqera.io/library/python_pip_pandas_plotly_upsetplot:1131451904bdd81f' }"

    input:
    tuple val(meta), path(input)

    output:
    path("*.html")              , emit: plots
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_svlendist.py \\
        $input \\
        -o ${params.variant_type}_lenght_${prefix}.mqc.html \\
        --tag ${meta.tag} \\
        --vartype ${params.variant_type} \\
        $args

    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${params.variant_type}_lenght_${prefix}.mqc.html

    """
}
