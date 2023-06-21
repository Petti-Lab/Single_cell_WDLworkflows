version 1.0

task run_clustering_n_pca_simple{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File clustering_script
    File rds_file_path
    String Sample_name
    Float clustering_res
    Int nPC
    File tirosh_file_path
  }

   command <<<
    Rscript ~{clustering_script} ~{rds_file_path} ~{Sample_name} ~{clustering_res} ~{nPC} ~{tirosh_file_path}
    >>>

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
  File clustering_rds = glob("*.Cycling.SCT.PCA.UMAP.TSNE.CLUST.*.rds")[0]
  }
}

workflow Seurat_clustering_simple {

  call run_clustering_n_pca_simple
}
