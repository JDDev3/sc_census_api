 
#devtools::install_github("cellgeni/schard")
# convert h5ad dataset to seurat
test_seurat <- schard::h5ad2seurat("test.h5ad")
# pseudobulk data to each donor
test_seurat_aggregated <- Seurat::AggregateExpression(test_seurat, return.seurat = TRUE, group.by = "donor_id")
# extract expression data
test_seurat_aggregated_data <- as.data.frame(test_seurat_aggregated@assays$RNA$data)
# extract gene information
test_seurat_aggregated_data$Ensemble <- test_seurat@assays$RNA@meta.features$feature_id
test_seurat_aggregated_data$Genes <- test_seurat@assays$RNA@meta.features$feature_name
test_seurat_aggregated_data <- test_seurat_aggregated_data |> dplyr::relocate("Ensemble")
test_seurat_aggregated_data <- test_seurat_aggregated_data |> dplyr::relocate("Genes")
rownames(test_seurat_aggregated_data) <- test_seurat_aggregated_data$Genes
# remove genes with all zero values across all samples
test_seurat_aggregated_data_nozero <- test_seurat_aggregated_data[apply(test_seurat_aggregated_data[,-c(1,2)], 1, function(x) !all(x == 0)),]
# extract meta information
test_seurat_aggregated_meta <- as.data.frame(test_seurat_aggregated@meta.data)
# write tsv files for expression and meta data 
readr::write_tsv(test_seurat_aggregated_data_nozero, "test_seurat_aggregated_data_nozero.tsv")
readr::write_tsv(test_seurat_aggregated_meta, "test_seurate_aggregated_meta.tsv")
