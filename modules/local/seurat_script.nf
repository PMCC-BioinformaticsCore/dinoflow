process SEURAT {
    label 'process_single'

    container "rlupat/r-dinoflow:0.1.1"
    
    input:
    tuple val(meta), path(anno), path(solo_dir)
    
    output:
    path "*.rds", emit: rds 
    path "*.csv", emit: annotation

    script:
    matrix   = "${solo_dir}/Gene*/filtered/matrix.mtx.gz"
    barcodes = "${solo_dir}/Gene*/filtered/barcodes.tsv.gz"
    features = "${solo_dir}/Gene*/filtered/features.tsv.gz"   
    output_folder = "QC_out/"
    out_RDS_file = "seuratObj.rds"
    out_ann_file = "annotation-mod.csv"
 
    """
    if [[ ! -d ${output_folder} ]]; then mkdir -p ${output_folder}; fi 

    seurat_script.R \\
        ${meta.id} \\
        ${anno} \\
        ${matrix} \\
        ${barcodes} \\
        ${features} \\
        ${output_folder} \\
        ${out_RDS_file} \\
        ${out_ann_file} 
    """
}
