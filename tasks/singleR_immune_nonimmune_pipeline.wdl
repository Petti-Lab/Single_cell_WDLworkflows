version 1.0

task run_singleR_immune_nonimmune {

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File singleR_immune_nonimmune_rscript
    File rds_file_path
    String Sample_name
    Float clustering_res
    Int nPC
    File singleR_yaml
  }

   command <<<
    Rscript ~{singleR_immune_nonimmune_rscript} ~{rds_file_path} ~{Sample_name} ~{clustering_res} ~{nPC} ~{singleR_yaml}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_rds = glob("*.rds")
  Array[File] output_jpgs = glob("*.jpg")
  }
}

workflow singleR_immune_nonimmune {

  call run_singleR_immune_nonimmune
}
