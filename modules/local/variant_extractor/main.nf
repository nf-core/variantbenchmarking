process VARIANT_EXTRACTOR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c9/c9964f2ad344f5df543a390ddceb308af0816bed9c9abb9cb400e65525580694/data':
        'community.wave.seqera.io/library/htslib_pysam_tabix_pip_variant-extractor:9a7098708edb08d4' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.norm.vcf.gz")   , emit: output
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    normalize_vcf.py \\
        $input \\
        ${prefix}.norm.vcf \\
        $fasta

    bgzip ${prefix}.norm.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.norm.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
