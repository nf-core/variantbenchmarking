process WITTYER {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://kubran/wittyer:0.3.3.0' :
        'kubran/wittyer:0.3.3.0' }"

    input:
    tuple val(meta),path(vcf), path(tbi), path(truth_vcf), path(truth_tbi), path(bed)
    path(config)

    output:
    tuple val(meta),    path("*ConfigFileUsed.json") , emit: config
    tuple val(meta),    path("*.Stats.json")         , emit: report
    tuple val(meta),    path("*eval.vcf.gz")         , emit: bench_vcf
    tuple val(meta),    path("*eval.vcf.gz.tbi")     , emit: bench_vcf_gzi
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = bed ? "--includeBed=$bed" : ""
    def config = config ? "--configFile=$config" : ""

    """
    mkdir bench
    dotnet /opt/Wittyer/Wittyer.dll \\
        --truthVcf=${truth_vcf} \\
        --inputVcf=${vcf} \\
        --outputDirectory=bench \\
        ${regions} \\
        ${config} \\
        ${args}

    mv bench/Wittyer.ConfigFileUsed.json ${prefix}.ConfigFileUsed.json
    mv bench/Wittyer.Stats.json ${prefix}.Stats.json
    mv bench/*.vcf.gz ${prefix}.eval.vcf.gz
    mv bench/*.vcf.gz.tbi ${prefix}.eval.vcf.gz.tbi

    rm -rf bench

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wittyer: 0.3.3.0
    END_VERSIONS
    """

}
