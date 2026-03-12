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
    tuple val("${task.process}"), val('cat'), eval("cat --version | head -n 1 | sed 's/^.*(GNU coreutils) //; s/ Copyright.*\$//'"), emit: versions_cat, topic: versions

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