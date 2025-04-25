process BEDTOOLS_INTERSECT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f12bf7e05df84523fa897da87c135edef2dee90a0567a6bf70f2f81ada6cc25/data' :
        'community.wave.seqera.io/library/pybedtools_pandas:6eae56b7891a94c6' }"

    input:
    tuple val(meta), path(truth), path(test)

    output:
    tuple val(meta),path("*stats.csv")    , emit: summary
    tuple val(meta),path("*TP.bed")       , emit: tp
    tuple val(meta),path("*FP.bed")       , emit: fp
    tuple val(meta),path("*FN.bed")       , emit: fn
    tuple val(meta),path("*converted.bed"), emit: out_bed, optional: true
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // caveman, cnvkit, controlfreec, facets and ascat has different cnv file beds which needs to be converted, yet if a
    def format = (meta.caller == "caveman" || meta.caller == "cnvkit" || meta.caller == "controlfreec" || meta.caller == "facets" || meta.caller == "ascat") && (!meta.converted) ? "${meta.caller}" : "bed"

    """
    bedtools_intersect.py \\
        $truth \\
        $test \\
        $prefix \\
        $format \\
        $params.genome

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
