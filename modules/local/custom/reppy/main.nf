process REPPY {
    tag "reppy"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/ebaae67bf8eaff1be0b2cde494309e8a8d0b7ae3778acbd54d47c10652e95990/data':
        'community.wave.seqera.io/library/jinja2:3.1.6--15077719b52e3f45' }"

    input:
    path (roc_csv)
    val  (reppy_input)

    output:
    path("*.html"), emit: html
    tuple val("${task.process}"), val('reppy'), val("0.1.0 "), topic: versions, emit: versions_reppy

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix 

    """
    rep.py \\
        $args \\
        $reppy_input \\
        -o ${prefix}.report.html \\

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix

    """
    echo $args
    
    touch ${prefix}.report.html
    """
}
