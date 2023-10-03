process EXTRACT_BARCODE_TXT {
    label 'process_single'

    conda 'ubuntu:20.04'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'biocontainers/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reads)

    output:
    path 'barcode.txt', emit: barcode

    script:
    """
    cut -f\$',' -d 9 ${meta.annotation} | tail -n+2 > barcode.txt
    """
}