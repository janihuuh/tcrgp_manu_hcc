
print("Create Seurat object...")

## Create Seurat object
counts           <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
hcc_seurat       <- CreateSeuratObject(counts, min.cells = 3,  min.features = 200, project = "hcc", meta.data = tcrab)



## Use only CD8 cells for this purpose
cd8_use         <- hcc_seurat@meta.data$Cell.Name[hcc_seurat@meta.data$phenotype == "CD8"]
hcc_cd8         <- subset(hcc_seurat, cells = cd8_use)
Idents(hcc_cd8) <- hcc_cd8@meta.data$T.cell.Cluster


## Log-normalise and scale the data and look for HVGs
hcc_cd8 <- NormalizeData(hcc_cd8, normalization.method = "LogNormalize", scale.factor = 10000)
hcc_cd8 <- FindVariableFeatures(object = hcc_cd8, mean.function = ExpMean, dispersion.function = LogVMR,
                             x.low.cutoff = 0.5, x.high.cutoff = 3, y.cutoff = 0.5,
                             do.plot = F)
hcc_cd8 <- ScaleData(hcc_cd8, do.par = TRUE, num.cores = 3, display.progress = T)



## Run linear dimensionality reduction
## Select only the PC:s with sdev > 2
hcc_cd8 <- RunPCA(object = hcc_cd8, pc.genes = VariableFeatures(hcc_cd8), pcs.compute = 100, do.print = FALSE)
nPC     <- sum(hcc_cd8[['pca']]@stdev > 2)



## Run non-linear dimensionality reduction
# hcc_cd8 <- RunUMAP(object = hcc_cd8, reduction.use = "pca", dims = 1:nPC, learning.rate = 1) ## Doesn't work for some reason, use a work-around
umap_df           <- umap::umap(hcc_cd8[['pca']]@cell.embeddings[,1:nPC])
umap_df           <- umap_df$layout %>% as.data.frame()
colnames(umap_df) <- c('UMAP1', 'UMAP2')

# hcc_cd8@meta.data <- cbind(hcc_cd8@meta.data, umap_df)
viz_df            <- hcc_cd8@meta.data
viz_df            <- cbind(viz_df, umap_df)
viz_df$cluster    <- viz_df$T.cell.Cluster

write.table(viz_df, paste0(output_dir, "viz_df.txt"), sep = "\t", quote = F, row.names = F)


## Save
saveRDS(hcc_cd8, paste0(output_dir, "hcc_cd8_seurat_object.rds"))
