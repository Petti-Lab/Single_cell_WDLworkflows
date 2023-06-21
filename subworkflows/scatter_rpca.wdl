version 1.0

import "../tasks/rpca_integration.wdl" as rpca


workflow scatter_rpca{

	input {
	String docker_image
    String queue_name
    Int mem_gb
    File seurat_rpca_rscript
    File inputSamplesFile
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	}

	scatter (sample in inputSamples) {
		call rpca.run_rpca_integration as run_rpca {
			input:
			    docker_image=docker_image,
    			queue_name=queue_name,
    			mem_gb=mem_gb,
    			seurat_rpca_rscript=seurat_rpca_rscript,
    			Seurat_rds=sample[1],
    			output_prefix=sample[0],
    			subset_column=sample[2],
				subset_column_selection=sample[3]
	}}
}
