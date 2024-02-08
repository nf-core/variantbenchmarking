process SURVIVOR_STATS {
    tag "$meta.id $meta2.caller"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/survivor:1.0.7--h9a82719_1':
        'quay.io/biocontainers/survivor:1.0.7--h9a82719_1' }"

    input:
    tuple val(meta),val(meta2), path(vcf), path(index)
    val(minsv)          // Min SV size (-1 to disable)
    val(maxsv)          // Max SV size (-1 to disable)
    val(minnumreads)    // Min number of reads support: RE flag (-1 to disable)

    output:
    tuple val(meta),val(meta2), path("*.stats"), emit: stats
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def name   = vcf.getBaseName()
    
    """
    gzip -d $vcf 

    SURVIVOR \\
        stats \\
        $name \\
        $minsv \\
        $maxsv \\
        $minnumreads \\
        ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.stats
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        survivor: \$(echo \$(SURVIVOR 2>&1 | grep "Version" | sed 's/^Version: //'))
    END_VERSIONS
    """
}
