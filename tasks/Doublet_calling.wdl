version 1.0

task run_doublet_collection{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File doublet_calling_script
    String cellranger_outs_directory
    File barcodes = "~{cellranger_outs_directory + '/'}barcodes.tsv.gz"
    File features = "~{cellranger_outs_directory + '/'}features.tsv.gz"
    File matrix = "~{cellranger_outs_directory + '/'}matrix.mtx.gz"
    String Sample_name
  }

   command <<<
    mkdir ~{Sample_name}
    mkdir ~{Sample_name}/filtered_feature_bc_matrix
    cp ~{barcodes} ~{Sample_name}/filtered_feature_bc_matrix/
    cp ~{features} ~{Sample_name}/filtered_feature_bc_matrix/
    cp ~{matrix} ~{Sample_name}/filtered_feature_bc_matrix/
    Rscript ~{doublet_calling_script} ~{Sample_name}/filtered_feature_bc_matrix/ ~{Sample_name}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File doublet_results = glob("*_Doublet_collection_results.rds")[0] 
  }
}

workflow Doublet_colletion{

  call run_doublet_collection
}

task run_seurat_singlesample{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File seurat_singlesample_rscript
    String cellranger_outs_directory
    File barcodes = "~{cellranger_outs_directory + '/'}barcodes.tsv.gz"
    File features = "~{cellranger_outs_directory + '/'}features.tsv.gz"
    File matrix = "~{cellranger_outs_directory + '/'}matrix.mtx.gz"
    String Sample_name
  }

   command <<<
    mkdir {Sample_name}/filtered_feature_bc_matrix
    cp ~{barcodes} {Sample_name}/filtered_feature_bc_matrix/
    cp ~{features} {Sample_name}/filtered_feature_bc_matrix/
    cp ~{matrix} {Sample_name}/filtered_feature_bc_matrix/
    Rscript ~{seurat_singlesample_rscript} {Sample_name}/filtered_feature_bc_matrix/ ~{Sample_name}