process BEDTOOLS_INTERSECT_BENCH {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/71/710294e18c04918cbd5dedd64abb3f43b34fb0f8399a1ed314316151a74a6147/data' :
        'community.wave.seqera.io/library/pybedtools_pandas:003b11d9e5690eba' }"

    input:
    tuple val(meta), path(truth), path(test)

    output:
    tuple val(meta),path("*stats.csv")    , emit: summary
    tuple val(meta),path("*TP_base.bed")  , emit: tp_base
    tuple val(meta),path("*TP_comp.bed")  , emit: tp_comp
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
        $params.genome \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    """
    touch ${meta.id}_stats.txt
    touch ${meta.id}_TP_comp.bed
    touch ${meta.id}_TP_base.bed
    touch ${meta.id}_FP.bed
    touch ${meta.id}_FN.bed
    touch ${meta.id}_converted.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
