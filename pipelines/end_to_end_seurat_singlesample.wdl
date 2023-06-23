version 1.0

import "../tasks/add_doubletinfo.wdl" as add_doubletinfo
import "../tasks/Doublet_calling.wdl" as Doublet_calling
import "../subworkflows/scatter_gather_singleR.wdl" as scatter_gather_singleR
import "../tasks/single_sample_seurat.wdl" as single_sample_filtering
import "../tasks/clustering_n_pca_simple.wdl" as single_sample_clustering
import "../subworkflows/recluster_renorm_rerun_singleR.wdl" as recluster_renorm_rerun_singleR

workflow end_to_end_seurat_single_sample{

    input {
    String docker_image
    String queue_name
    Int mem_gb
    File inputSamplesFile
    File singleR_refsFile
    String subset_column_name
    String ident_name
    String inverse
    String output_suffix
    String organism    
    File recluster_renormalize_script
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
    }

    scatter (sample in inputSamples) {
        call Doublet_calling.run_doublet_collection as run_doublet {
        input:
        docker_image=docker_image,
        queue_name=queue_name,
        mem_gb=mem_gb,
        Directory cellranger_outs_directory=sample[1],
        Sample_name=sample[0]
        }
        call single_sample_filtering.run_seurat_singlesample as run_single_srt {
        input:
        docker_image=docker_image,
        queue_name=queue_name,
        mem_gb=mem_gb,
        Directory cellranger_outs_directory=sample[1],
        Sample_name=sample[0]
        }
        call single_sample_clustering.run_clustering_n_pca_simple as run_srt_clust_simp {
        input:
        docker_image=docker_image,
        queue_name=queue_name,
        mem_gb=mem_gb,
        rds_file_path=run_single_srt.intermed_rds,
        Sample_name=sample[0]
        }
        call add_doubletinfo.add_doublets_metadata as add_doub_to_srt {
        input:
        docker_image=docker_image,
        queue_name=queue_name,
        mem_gb=mem_gb,
        input_rds_file=run_srt_clust_simp.clustering_rds,
        doublet_file=run_doublet.doublet_results
        }
        call scatter_gather_singleR.scatter_gather_singleR as scat_gath_singleR {
        input:
        docker_image=docker_image,
        queue_name=queue_name,
        mem_gb=mem_gb,
        seurat_rds=add_doub_to_srt.seurat_doublet_rds,
        inputSamplesFile=singleR_refsFile
        }
        call recluster_renorm_rerun_singleR.LinearChain_recluster_rerun_singleR as LinearChain_recluster_rerun_singleR {
        input:
        docker_image=docker_image,
        queue_name=queue_name,
        mem_gb=mem_gb,
        seurat_rds=scat_gath_singleR.seurat_singleR_rds,
        subset_column_name=subset_column_name,
        ident_name=ident_name,
        inverse=inverse,
        organism=organism,
        output_suffix=output_suffix,
        inputSamplesFile=singleR_refsFile
        }
}

 call merge_singlesample_seurat_to_multisample {
      input: 
        docker_image=docker_image, 
        queue_name=queue_name,
        mem_gb=mem_gb,
        seurat_files=LinearChain_recluster_rerun_singleR.seurat_singleR_rds,
    }
}


task merge_singlesample_seurat_to_multisample{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File merge_singlesample_seurat_script
    String organism
    String project_name
    Array[File] seurat_files
  }

   command <<<
    Rscript ~{merge_singlesample_seurat_script} ~{project_name} ~{organism} ~{sep=" " seurat_files}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File final_seurat_obj = glob("*.merged_multisample.*.rds")[0] 
  }
}
