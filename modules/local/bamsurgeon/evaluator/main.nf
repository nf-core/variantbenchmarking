process BAMSURGEON_EVALUATOR {
    tag "$meta.id"
    label 'process_low'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lethalfang/bamsurgeon:1.2':
        'lethalfang/bamsurgeon:1.2' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta),path("*falsepositives.vcf") , emit: fp
    tuple val(meta),path("*truepositives.vcf")  , emit: tp
    tuple val(meta),path("*falsenegatives.vcf") , emit: fn
    tuple val(meta),path("*.stats")             , emit: stats
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def muttype = meta.vartype.contains("snv") ? "SNV" : meta.vartype.contains("indel") ? "INDEL" : meta.vartype.contains("sv") ? "SV" : ""

    """
    python3 /usr/local/bamsurgeon/scripts/evaluator.py \\
        -v $vcf \\
        -t $truth_vcf \\
        -f $fasta \\
        -m $muttype \\
        $args \\
        --fp ${prefix}.falsepositives.vcf \\
        --tp ${prefix}.truepositives.vcf \\
        --fn ${prefix}.falsenegatives.vcf > ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamsurgeon: v1.2
    END_VERSIONS
    """
}
