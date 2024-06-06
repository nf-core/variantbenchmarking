process VISUALIZE_TAGGED {
    tag "$meta.benchmark_tool"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-shinyngs:1.8.4--r43hdfd78af_0':
        'biocontainers/r-shinyngs:1.8.4--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta),path("*.txt"), emit: plots
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"

    """
    grep -o 'SUPP_VEC=[^,;]*' $vcf | awk -F'=' '{print \$2}' | sed -e 's/\\(.\\)/\\1 /g' > ${vcf}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.benchmark_tool}"
    """
    touch ${vcf}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
