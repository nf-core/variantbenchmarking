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
    # Use Groovy to fill in the template and generate the YAML file
    # groovy -e \"
    def template = new File('\$template').text
    def engine = new groovy.text.SimpleTemplateEngine()
    def binding = [csvpath: '\$csv']
    def yaml = engine.createTemplate(template).make(binding).toString()
    new File('config.yaml').write(yaml)
    # \"
    """
}
