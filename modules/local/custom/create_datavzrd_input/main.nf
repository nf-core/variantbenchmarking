process CREATE_DATAVZRD_INPUT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch config.yaml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//')
    END_VERSIONS
    """
}
