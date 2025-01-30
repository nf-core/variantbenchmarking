process CREATE_DATAVZRD_INPUT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::tar=1.34"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(csv), path(template)

    output:
    tuple val(meta), path("*.yaml"), path(csv), emit: config

    script:
    """
    #!/bin/bash

    cat "$template" | sed "s|CSVPATH|$csv|g" > config.yaml
    """

    stub:
    """
    touch config.yaml
    """
}
