process EXTRACT_MAIN {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.txt"), emit: chr_list
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def name = bed.getBaseName()

    """
    cat $bed | awk '{print \$1}' | uniq > ${name}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
