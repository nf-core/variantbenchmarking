process ONCOLINER_ASSESSMENT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pysam_jinja2_libzio_libzip_pruned:4032adb7b5c672a5':
        'community.wave.seqera.io/library/pysam_jinja2_libzio_libzip_pruned:4032adb7b5c672a5' }"

    input:
    tuple val(meta), path(test), path(truth), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*")   , emit: output
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = bed ? "--bed-masks $bed" : ""

    """
    assessment_main.py \\
        -t $truth \\
        -v $test \\
        -f $fasta \\
        -o ${prefix} \\
        $bed \\
        --keep-intermediates

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
