process VCF_GENOTYPE_ANNOTATOR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://griffithlab/vatools:5.1.10'
        'griffithlab/vatools:5.1.10' }"

    input:
    tuple val(meta),path(vcf)

    output:
    tuple val(meta),path("*.{vcf}"), emit: vcf
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    """
    vcf-genotype-annotator \\
        $vcf \\
        ${meta.id} \\
        $args
        -o ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
