process SEURAT_TO_EDGER {
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
    path "*smearPlots_out.gz"
    
    """
    Rscript seurat_script.R \\
        ${seurat_obj_file_path} \\
        ${mod_annotation_file_path} 
    gzip smearPlots_out/
    """
}
