process ADDHEAD {
    tag "$meta.id $meta2.caller" 
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.18--h8b25389_0':
        'quay.io/biocontainers/bcftools:1.18--h8b25389_0' }"

    input:
    tuple val(meta), val(meta2), path(vcf), path(header)

    output:
    tuple val(meta), val(meta2), path("*.vcf.gz")  , emit: vcf
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cat $header $vcf > ${prefix}_rehead.vcf

    bcftools \\
        sort \\
        --output ${prefix}.vcf.gz \\
        --output-type z \\
        --temp-dir . \\
        ${prefix}_rehead.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
