version 1.0

task run_scsorter {

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File scsorter_script
    String rds_file_path
    String Sample_name
    String marker_list_rds
    String marker_list_name
  }

   command <<<
    Rscript ~{scsorter_script} ~{rds_file_path} ~{Sample_name} ~{marker_list_rds} ~{marker_list_name}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_rds = glob("*.rds")
  }
}

workflow scsorter_markerprediction {

  call run_scsorter
}