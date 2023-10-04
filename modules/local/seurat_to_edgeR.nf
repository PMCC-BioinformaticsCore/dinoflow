process MTX_TO_SEURAT {
    //tag "$meta.id"
    //label 'process_medium'

    container "rlupat/r-dinoflow:0.1.1"
    
    publishDir params.outdir
    
    input:
  
    path seurat_obj_file_path
    path mod_annotation_file_path
    
    output:
    
    path "*BCV_plot.png"
    path "*CSV_out/*_table.csv"
    path "*smearPlots_out.zip"
    
    """
    Rscript seurat_script.R \\
        ${seurat_obj_file_path} \\
        ${mod_annotation_file_path} 
    gzip ${output_plot_smear_dir}
    """
}