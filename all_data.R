suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(SingleCellExperiment)))

output_path    <- "final_matrices/"

all_files      <- list.files("../MarkerPaper/QCfiltered_data")
data_files     <- all_files[grepl(".rds", all_files, fixed=TRUE)]

common_genes   <- readRDS("../MarkerPaper/rds_files/common_genes_agg.rds")

final_matrix   <- c()
final_metadata <- c()

for (data in data_files){
	
	print(data)
	
	matrix             <- readRDS(paste0("../MarkerPaper/sces/", data))

	counts             <- counts(matrix)
	metadata           <- colData(matrix)

	colnames(counts)   <- paste0(colnames(counts),   "_", tstrsplit(data, ".rds", fixed=TRUE)[[1]])
	rownames(metadata) <- paste0(rownames(metadata), "_", tstrsplit(data, ".rds", fixed=TRUE)[[1]]) 

	final_matrix       <- cbind(final_matrix, counts[common_genes, ])

	df           <- as.data.frame(metadata[, "cluster_id"])
	rownames(df) <- rownames(metadata)
	colnames(df) <- "Cell"

	final_metadata     <- rbind(final_metadata, df)

}

saveRDS(final_matrix,   paste0(output_path, "all_data.rds"))
saveRDS(final_metadata, paste0(output_path, "all_metadata.rds"))
