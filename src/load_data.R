
set.seed(123)

library(ggplot2, verbose = F)
library(dplyr, verbose = F)
library(data.table, verbose = F)
library(RColorBrewer, verbose = F)
library(Seurat, verbose = F)

add_guide   <- guides(colour = guide_legend(override.aes = list(size=5)))

getPalette  <- colorRampPalette(brewer.pal(7, "Set2"))
getPalette2 <- colorRampPalette(brewer.pal(5, "Set1"))

## Read in the TCRab file
tcrab  <- fread("data/zheng_tcrab.txt")

## Read in the data. We start with counts as they can be the most reliably modelled
# Count data downloaded from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98638
counts    <- fread("data/zheng_count.txt") %>% as.data.frame()
cellnames <- colnames(counts)
genenames <- rownames(counts)

## Read in the tcrgp-data
tcrgp  <- fread("data/TCRGP/HCC_b_85.csv")
