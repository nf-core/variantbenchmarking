process VCF_TO_CSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9964f2ad344f5df543a390ddceb308af0816bed9c9abb9cb400e65525580694/data':
        'community.wave.seqera.io/library/htslib_pysam_tabix_pip_variant-extractor:a12ef217eccf6ba8' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.csv")   , emit: output
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vcf_to_csv.py \\
        $input \\
        ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
