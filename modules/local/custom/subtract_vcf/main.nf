process SUBTRACT_VCF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/ab3b0054e3111812d8f2deb12345d5b7ca7ea7b18a2dbcbf174d46274c28deba/data':
        'community.wave.seqera.io/library/pip_pandas:40d2e76c16c136f0' }"

    input:
    tuple val(meta), path(vcf), path(index), path(regions), path(targets), path(exclude)

    output:
    tuple val(meta), path("*remain.vcf.gz"), emit: vcf, optional:true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = regions ? "--bed-file ${regions}" : ''

    """
    subtract_variants.py \\
        $vcf \\
        $exclude \\
        ${prefix}.remain.vcf.gz \\
        $bed \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.remain.vcf.gz
    """
}