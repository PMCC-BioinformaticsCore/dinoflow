process MTX_TO_SEURAT {
    //tag "$meta.id"
    //label 'process_medium'

    conda "r-seurat"
    container "nf-core/seurat:4.3.0"
    
    publishDir params.outdir
    
    input:
    val prefix
    path ann_data
    path matrix
    path barcodes
    path features
    path output_folder
    path out_RDS_file
    path out_ann_file
    
    output:
    path $output_folder
    
    """
    Rscript seurat_script.R \\
        ${prefix} \\
        ${ann_data} \\
        ${matrix} \\
        ${barcodes} \\
        ${features} \\
        ${output_folder} \\
        ${out_RDS_file} \\
        ${out_ann_file} 
    """
}
