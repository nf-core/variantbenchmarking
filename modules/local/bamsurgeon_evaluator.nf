process BAMSURGEON_EVALUATOR {
    tag "$meta.id"
    label 'process_low'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lethalfang/bamsurgeon:1.2':
        'lethalfang/bamsurgeon:1.2' }"

    input:
    tuple val(meta),val(meta2), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi)
    tuple path(fasta), path(fai)
    val(muttype)

    output:
    tuple val(meta),val(meta2), path("*.{vcf}"), emit: bench
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    python3 /usr/local/bamsurgeon/scripts/evaluator.py \\
        -v $vcf \\
        -t $truth_vcf \\
        -f $fasta \\
        -m $muttype \\
        $args \\
        --fp ${prefix}.falsepositives.vcf \\
        --tp ${prefix}.truepositives.vcf \\
        --fn ${prefix}.falsenegatives.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamsurgeon: v1.2
    END_VERSIONS
    """
}
