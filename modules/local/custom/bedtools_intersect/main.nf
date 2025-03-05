process BEDTOOLS_INTERSECT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ab/ab3b0054e3111812d8f2deb12345d5b7ca7ea7b18a2dbcbf174d46274c28deba/data':
        'community.wave.seqera.io/library/pip_pandas:40d2e76c16c136f0' }"

    input:
    tuple val(meta), path(truth), path(test)

    output:
    tuple val(meta),path("*stats.txt")    , emit: summary
    tuple val(meta),path("*TP.bed")       , emit: tp
    tuple val(meta),path("*FP.bed")       , emit: fp
    tuple val(meta),path("*FN.bed")       , emit: fn
    tuple val(meta),path("*converted.bed"), emit: out_bed, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def format = meta.caller == "caveman" || meta.caller == "cnvkit" || meta.caller == "controlfreec" || meta.caller == "facets" ? "${meta.caller}" : "bed"

    """
    bedtools_intersect.py
        $truth
        $test
        $meta.id
        $format

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}_stats.txt
    touch ${meta.id}_TP.bed
    touch ${meta.id}_FP.bed
    touch ${meta.id}_FN.bed
    touch ${meta.id}_converted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
