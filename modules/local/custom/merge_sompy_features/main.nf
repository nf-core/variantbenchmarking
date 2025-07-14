process MERGE_SOMPY_FEATURES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/ab3b0054e3111812d8f2deb12345d5b7ca7ea7b18a2dbcbf174d46274c28deba/data':
        'community.wave.seqera.io/library/pip_pandas:40d2e76c16c136f0' }"

    input:
    tuple val(meta), path(csvs)

    output:
    tuple val(meta), path("*.summary.csv")   , emit: TP
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    merge_sompy_features.py $csvs --output ${prefix}.${meta.tag}.summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${meta.tag}.summary.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
