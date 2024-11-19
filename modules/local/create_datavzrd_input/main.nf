process CREATE_DATAVZRD_INPUT {
    tag "$meta.id"
    label 'process_single'

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
