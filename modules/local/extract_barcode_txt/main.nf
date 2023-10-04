process EXTRACT_BARCODE_TXT {
    label 'process_single'

    conda 'ubuntu:20.04'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/ubuntu:20.04' }"

    input:
    path(anno)

    output:
    path 'barcode.txt', emit: barcode

    shell:
    """ 
    rev !{anno} | awk -F "," '{print \$1}' | tr -d '"' | tail -n+2 > barcode.txt
    """
}
