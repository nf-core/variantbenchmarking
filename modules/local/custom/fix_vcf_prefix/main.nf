process FIX_VCF_PREFIX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d097f8fee4db1239705e5d3929e50547796ab91a1e0bdf4201030ced8cb272d/data':
        'community.wave.seqera.io/library/bcftools_pip:8a460830487271c2' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(rename_chr)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"

    if ("$input" == "${prefix}.vcf.gz") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    fix_vcf_prefix.py \\
        $input \\
        ${prefix}.vcf.gz \\
        --rename-chr $rename_chr \\
        --target-version $params.genome

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo '' | gzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
