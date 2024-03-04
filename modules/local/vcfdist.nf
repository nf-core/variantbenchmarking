process VCFDIST {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://timd1/vcfdist:v2.3.2' :
        'timd1/vcfdist:v2.3.2' }"

    input:
    tuple val(meta),path(vcf), path(tbi), path(truth_vcf), path(truth_tbi), path(bed)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta), path("*.tsv,vcf"), emit: bench
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "-b $bed" : ""

    """
    vcfdist \\
        ${vcf} \\
        ${truth_vcf} \\
        $fasta \\
        -p ${prefix} \\
        ${regions} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcfdist: \$(echo \$(vcfdist --version 2>&1) | sed 's/^.*vcfdist v//')
    END_VERSIONS
    """

}
