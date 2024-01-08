process EXTRACT_BARCODE_TXT {
    label 'process_single'

    container "rlupat/r-dinoflow:0.1.1"

    input:
    path(anno)

    output:
    path 'barcode.txt', emit: barcode

    shell:
    """ 
    Rscript -e 'data = read.csv("$anno",sep=","); cat(data[,"Barcode"])' > barcode.txt
    """
}
