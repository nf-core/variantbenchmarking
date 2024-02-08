process WITTYER {
    tag "$meta.id $meta2.caller"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/wittyer:0.3.3.0' :
        'kubran/wittyer:0.3.3.0' }"

    input:
    tuple val(meta),val(meta2), path(vcf), path(tbi), path(truth_vcf), path(truth_tbi), path(bed)
    path(config)

    output:
    tuple val(meta),    path(wittyer_bench) , emit: bench
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--includebed $bed" : ""
    def config = bed ? "--configFile $config" : ""

    """
    dotnet Wittyer.dll \\
        --truthVcf ${truth_vcf} \\
        --inputVcf ${vcf} \\
        --output wittyer_bench \\
        ${regions} \\
        ${config} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wittyer: \$(echo \$(dotnet Wittyer.dll --version 2>&1) | sed 's/^.*witty\.er v//')
    END_VERSIONS
    """
    
}
