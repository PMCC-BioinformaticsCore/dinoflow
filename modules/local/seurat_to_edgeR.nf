process SEURAT_TO_EDGER {
    container "rlupat/r-dinoflow:0.1.1"
    
    input:
  
    path seurat_obj_file_path
    path mod_annotation_file_path
    
    output:
    
    path "*BCV_plot.png"
    path "*CSV_out/*_table.csv"
    path "*smearPlots_out.gz"
    
    """
    seurat_script.R \\
        ${seurat_obj_file_path} \\
        ${mod_annotation_file_path} 
    gzip smearPlots_out/
    """
}
