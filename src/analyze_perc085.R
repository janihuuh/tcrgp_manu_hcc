
## Analyze with 0.85 percent prob
global_folder <- here::here()
setwd(global_folder)

## The output folder
output_dir <- "results/"
dir.create(output_dir, showWarnings = F)

## ====== Run the analyses

## Initialize; download the data
source("src/fun_functions.R")
source("src/run_installPackages.R")
source("src/load_data.R")
source("src/run_modifyTCRGP.R")

## Create Seurat-object; do clustering
source("src/run_createSeurat.R")
source("src/run_clustering.R")

## Choose the new clustering
Idents(hcc_cd8) <- plyr::revalue(as.factor(hcc_cd8@meta.data$RNA_snn_res.0.6),  replace = c("0" = "exhausted 1", "1" = "memory", "2" = "effector", "3" = "naïve", "4" = "exhausted 2", "5" = "exhausted 3"))
viz_df$cluster  <- plyr::revalue(as.factor(viz_df$RNA_snn_res.0.6), replace = c("0" = "exhausted 1", "1" = "memory", "2" = "effector", "3" = "naïve", "4" = "exhausted 2", "5" = "exhausted 3"))

## Analyze the results
source("src/run_fisher.R")
source("src/run_clonality.R")

## Visualize the results
source("src/run_plotManuscript.R")
