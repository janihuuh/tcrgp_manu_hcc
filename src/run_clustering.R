
message("Do clustering...")

## Do clustering
res     <- seq(0.1,3.1,0.5) %>% as.list()
nPC     <- sum(hcc_cd8[['pca']]@stdev > 2)

hcc_cd8 <- FindNeighbors(hcc_cd8, dims = 1:nPC)
hcc_cd8 <- FindClusters(hcc_cd8, resolution = res, random.seed = 123)

saveRDS(hcc_cd8, paste0(output_dir, "hcc_cd8.rds"))

## Add clusters to viz df
viz_df  <- cbind(hcc_cd8@meta.data, umap_df)
colnames(viz_df)  <- make.names(colnames(viz_df), unique=TRUE)


## Visualize
clusters <- paste0("RNA_snn_res.", res)

pdf(paste0(output_dir, "umap_different_seurat_clusters.pdf"), width = 8, height = 6)
lapply(clusters, FUN = function(x) {plotUMAP(viz_df = viz_df, cluster = viz_df[,x]) %>% print} )
dev.off()
