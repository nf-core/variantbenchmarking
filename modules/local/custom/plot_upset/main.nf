process PLOT_UPSET {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bcc3d7359ce8c6c53e98534593b4a8a91fec5c1ab4bccd66f39f11128c39a1c3/data' :
        'community.wave.seqera.io/library/pip_upsetplot_matplot_pandas:d9e1259bc972b7a4' }"

    input:
    tuple val(meta), path(files)

    output:
    path("*.png")         , emit: plot, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_upset.py \\
        --fp ${meta.id}.FP.csv \\
        --fn ${meta.id}.FN.csv \\
        --tp-base ${meta.id}.TP_base.csv \\
        --tp-comp ${meta.id}.TP_comp.csv \\
        --output ${prefix} \\
        --title "Upset plot for ${meta.id}"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_tp_fn_mqc.png
    touch ${prefix}_tp_fp_mqc.png
    """
}