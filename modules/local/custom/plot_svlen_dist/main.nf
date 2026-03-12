process PLOT_SVLEN_DIST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/24/24f902548e45e009de670f96def94e83a1da47af87e793389091413a0182a820/data' :
        'community.wave.seqera.io/library/matplotlib_numpy_pandas:1503a72c3e08341d' }"

    input:
    tuple val(meta), path(input)

    output:
    path("*.png")         , emit: plot
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_svlendist.py \\
        $input \\
        -o ${prefix}.${params.variant_type}.mqc.png \\
        --title "INDEL Length Distributions of ${meta.tag} Variants" \\
        $args
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.svlen.png
    """
}