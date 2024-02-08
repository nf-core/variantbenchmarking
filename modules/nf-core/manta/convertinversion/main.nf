process MANTA_CONVERTINVERSION {
    tag "$meta.id"
    label 'process_low'
    label 'error_retry'

    conda "bioconda::manta=1.6.0 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-40295ae41112676b05b649e513fe7000675e9b84:a0332aa38645fbb8969567731ce68cfb7f830ec4-0':
        'quay.io/biocontainers/mulled-v2-40295ae41112676b05b649e513fe7000675e9b84:a0332aa38645fbb8969567731ce68cfb7f830ec4-0' }"

    input:
    tuple val(meta),val(meta2), path(vcf), path(index)
    tuple path(fasta), path(fai)

    output:
    tuple val(meta),val(meta2), path("*.vcf.gz"),path("*.vcf.gz.tbi")  , emit: vcf_tabi
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta2.caller == "manta"){
        """
        convertInversion.py \$(which samtools) $fasta $vcf | bgzip --threads $task.cpus > ${prefix}.converted.vcf.gz
        tabix ${prefix}.converted.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            manta: \$( configManta.py --version )
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
        END_VERSIONS
        """
    }
    else{
        """
        cp $vcf ${prefix}.vcf.gz
        tabix ${prefix}.vcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            manta: \$( configManta.py --version )
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
        END_VERSIONS
        """       
    }

}
