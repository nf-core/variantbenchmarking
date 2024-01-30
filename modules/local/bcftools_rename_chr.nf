process BCFTOOLS_RENAME_CHR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0':
        'quay.io/biocontainers/bcftools:1.17--haef29d1_0' }"

    input:
    tuple val(meta),val(meta2), path(input), path(index)
    each path(rename_file)

    output:
    tuple val(meta),val(meta2), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    "vcf"
    if ("$input" == "${prefix}.${extension}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    bcftools \\
        annotate \\
        $args \\
        --rename-chrs $rename_file \\
        --output ${prefix}.${extension} \\
        --threads $task.cpus \\
        --force \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
