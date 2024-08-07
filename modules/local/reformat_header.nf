process REFORMAT_HEADER {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(index)

    output:
    tuple val(meta),path("*reformatted.vcf.gz"), path("*reformatted.vcf.gz.tbi"), emit: gz_tbi
    path  "versions.yml"                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zcat $vcf | \
        awk 'BEGIN {OFS="\t"} /^##FORMAT=<ID=PS,/ {sub(/Type=Integer/, "Type=String")} {print}' | \
        sed 's/##FORMAT=<ID=AD,Number=[^,]*/##FORMAT=<ID=AD,Number=./' | \
        bgzip -c > ${prefix}.reformatted.vcf.gz

    tabix -p vcf ${prefix}.reformatted.vcf.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reformatted.vcf.gz
    touch ${prefix}.reformatted.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
