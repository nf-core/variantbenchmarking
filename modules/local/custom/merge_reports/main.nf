process MERGE_REPORTS {
    tag "$meta.benchmark_tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/ab3b0054e3111812d8f2deb12345d5b7ca7ea7b18a2dbcbf174d46274c28deba/data':
        'community.wave.seqera.io/library/pip_pandas:40d2e76c16c136f0' }"

    input:
    tuple val(meta), path(inputs)

    output:
    tuple val(meta),path("*.summary.csv")    , emit: summary
    tuple val(meta),path("*.regions.csv")    , emit: regions, optional: true
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //g'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"
    """
    merge_reports.py $inputs \\
        -b $meta.benchmark_tool \\
        -v $meta.vartype \\
        -a $params.analysis \\
        -o ${prefix}

    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"
    """
    touch ${prefix}.summary.csv
    touch ${prefix}.regions.csv

    """

}
