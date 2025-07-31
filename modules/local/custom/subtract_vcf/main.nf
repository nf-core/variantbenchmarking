process SUBTRACT_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1c/1cc0a8019c9e5473261c959c2ca234914fa61267cd4bd7c66a2f8cbfc7e06f70/data' :
        'community.wave.seqera.io/library/bcftools_bedtools:8ae45183a3924432' }"

    input:
    tuple val(meta), path(vcf), path(index), path(regions), path(targets), path(exclude)

    output:
    tuple val(meta), path("*remain.vcf.gz"), emit: vcf, optional:true
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env bash
    zcat $exclude | awk 'BEGIN{OFS="\\t"} {print \$1, \$2-1, \$2}' > exclude.bed

    if [ -s "${regions}" ]; then
        bcftools view -R ${regions} $vcf -O z -o ${prefix}.truth.vcf.gz
    elif [ -s "${targets}" ]; then
        bcftools view -T ${targets} $vcf -O z -o ${prefix}.truth.vcf.gz
    else
        mv $vcf ${prefix}.truth.vcf.gz
    fi

    if [ -s exclude.bed ]; then
        bedtools intersect -v -a ${prefix}.truth.vcf.gz -b exclude.bed > ${prefix}.remain.vcf
        gzip -c ${prefix}.remain.vcf > ${prefix}.remain.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
