process CREATE_DATAVZRD_INPUT {
    tag "$meta.id"
    label 'process_single'

    input:
    path(template) //tuple val(meta),
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.yaml"), path(csv), emit: config

    script:
    """
    #!/bin/bash

    cat "$template" | sed "s|CSVPATH|$csv|g" > config.yaml
    """
}
