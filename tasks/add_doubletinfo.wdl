version 1.0

task add_doublets_metadata{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File add_meta_script
    File doublet_file
    File input_rds_file
  }

   command <<<
    Rscript ~{add_meta_script} ~{input_rds_file} ~{doublet_file}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File seurat_doublet_rds = glob("*.doublet_calls.*rds")[0] 
  }
}

workflow Add_doublet_info{

  call add_doublets_metadata
}
