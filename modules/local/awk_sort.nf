process AWK_SORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta),val(meta2), path(vcf)

    output:
    tuple val(meta),val(meta2), path("*.tmp.vcf.gz"),path("*.tmp.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def zipname  = vcf.getBaseName()

    if (vcf.getExtension() == "gz"){
    """
    bgzip -d $vcf
    cat $zipname | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > ${zipname}.tmp.vcf
    bgzip -c ${zipname}.tmp.vcf > ${zipname}.tmp.vcf.gz
    tabix -p vcf ${zipname}.tmp.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
    }
    else{
    """
    cat $vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > ${zipname}.tmp.vcf
    bgzip -c ${zipname}.tmp.vcf > ${zipname}.tmp.vcf.gz
    tabix -p vcf ${zipname}.tmp.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
    }
}
