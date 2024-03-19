// UNTESTED
process GRIDSS_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gridss:2.13.2--h270b39a_0':
        'biocontainers/gridss:2.13.2--h270b39a_0' }"

    input:
    tuple val(meta), path(vcf), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta),path("*.vcf.gz"),path("*.vcf.gz.tbi")   , emit: vcf
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def genome = params.genome.contains("38") ? "hg38": "hg19"
    def VERSION = '2.13.2' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.

    """
    bgzip -d $vcf -c > unzziped.vcf
    simple_event-annotator.R \\
        unzziped.vcf \\
        ${prefix}.vcf \\
        ${genome}

    bgzip --threads ${task.cpus} -c ${prefix}.vcf > ${prefix}.anno.vcf.gz
    tabix -p vcf ${prefix}.anno.vcf.gz

    rm unzziped.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: ${VERSION}
    END_VERSIONS
    """

}
