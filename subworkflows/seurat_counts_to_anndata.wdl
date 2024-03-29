version 1.0

task Get_seurat_counts{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File Get_counts_script
    File seurat_rds
    String Sample_name
  }

   command <<<
    Rscript ~{Get_counts_script} ~{seurat_rds} ~{Sample_name}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File counts_matrix = glob("*_RNA_assay_counts_matrix_transposed*.txt")[0] 
  }
}

task make_anndata_object{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File make_anndata_py
    File counts_matrix_file
    String Sample_name
    String Organism
  }

   command <<<
    export NUMBA_CACHE_DIR="$PWD"
    python3 ~{make_anndata_py} ~{counts_matrix_file} ~{Sample_name} ~{Organism}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File anndata_object = glob("anndata_*.h5ad")[0] 
  }
}

workflow seurat_counts_to_anndata{

  input {
    String queue_name
    Int mem_gb
    String Sample_name
  }

  call Get_seurat_counts {
    input:
    queue_name=queue_name,
    mem_gb=mem_gb,
    Sample_name=Sample_name
  }
  call make_anndata_object {
    input:
    queue_name=queue_name,
    mem_gb=mem_gb,
    Sample_name=Sample_name,
    counts_matrix_file=Get_seurat_counts.counts_matrix,
  }
}