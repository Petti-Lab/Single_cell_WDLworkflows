version 1.0

task run_rpca_integration{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File seurat_rpca_rscript
    String Seurat_rds
    String output_prefix
    String subset_column
    String subset_column_selection
  }

   command <<<
    Rscript ~{seurat_rpca_rscript} ~{Seurat_rds} ~{output_prefix} ~{subset_column} ~{subset_column_selection}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_jpgs = glob("*.pdf")
  Array[File] output_rds = glob("*.rds")
  Array[File] output_txts = glob("*.txt")
  }
}

workflow Seurat_rPCA_integration {

  call run_rpca_integration
}