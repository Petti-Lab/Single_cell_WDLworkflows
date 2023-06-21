version 1.0

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
      mkdir ~{Sample_name}
      mkdir ~{Sample_name}/filtered_feature_bc_matrix
      cp ~{barcodes} ~{Sample_name}/filtered_feature_bc_matrix/
      cp ~{features} ~{Sample_name}/filtered_feature_bc_matrix/
      cp ~{matrix} ~{Sample_name}/filtered_feature_bc_matrix/
      Rscript ~{seurat_singlesample_rscript} ~{Sample_name}/filtered_feature_bc_matrix/ ~{Sample_name}
  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_pdfs = glob("*.pdf")
  Array[File] output_pngs = glob("*.png")
  Array[File] output_rds = glob("*.rds")
  Array[File] output_jpgs = glob("*.jpg")
  File intermed_rds = glob("*.MergedFilteredSeuratObject.*.rds")[0] 
  }
}

workflow Seurat_single_sample {

  call run_seurat_singlesample
}

