


CalculateFisherEpitope  <- function(viz_df, epitope_temp, patient_temp, tissue_temp, stage_temp, alternative = "greater"){

  ## Count to which cluster the given TCRs are enriched to

  # @ params:
  # viz_df       = df, that contains the TCRGP predictions and variables to count enrichment on
  # epitope_temp = char, epitope to count enrichment to
  # patient_temp = char, patient to count enrichment to
  # tissue_temp  = char, tissue to count enrichment to
  # alternative  = char, alternative for Fisher's exact test, 'less', 'greater' or 'two-sided'

  tot_df <- NULL
  i <- 1

  viz_df <- viz_df %>% filter(Patient %in% patient_temp) %>% filter(Sample.Type %in% tissue_temp) %>% filter(Stage %in% stage_temp)

  if(nrow(viz_df) == 0){
    # print('No samples')
    return()
  }

  patient_names <- paste(unique(viz_df$Patient), collapse=", ")
  tissue_names  <- paste(unique(viz_df$Sample.Type), collapse=", ")
  stage_names   <- paste(unique(viz_df$Stage), collapse=", ")

  epitope_names <- paste(epitope_temp, collapse=", ")

  ## Select cluster in the given epitope
  for(cluster_temp in unique(viz_df$cluster)){

    # print(cluster_temp)

    pred_df             <- viz_df %>% filter(cluster == cluster_temp) %>% select(epitope_temp)

    cluster_predicted   <- table(pred_df[,1])[2]
    cluster_unpredicted <- table(pred_df[,1])[1]

    epitope_df          <- viz_df %>% select(epitope_temp)

    predicted_cells     <- table(epitope_df[,1])[2] - cluster_predicted
    unpredicted_cells   <- table(epitope_df[,1])[1] - cluster_unpredicted

    if(is.na(cluster_predicted))   cluster_predicted = 0
    if(is.na(predicted_cells))     predicted_cells = 0
    if(is.na(cluster_unpredicted)) cluster_unpredicted = 0
    if(is.na(unpredicted_cells))   unpredicted_cells = 0

    df <- matrix(c(cluster_predicted, predicted_cells, cluster_unpredicted, unpredicted_cells), nrow = 2, byrow = T)

    fisher_result <- df %>%
      fisher.test(alternative = alternative) %>%
      broom::tidy()

    tot_df[[i]] <- cbind(patient_names,
                         epitope_names,
                         tissue_names,
                         stage_names,
                         cluster_temp,
                         cluster_predicted,
                         predicted_cells,
                         cluster_unpredicted,
                         unpredicted_cells,
                         fisher_result)
    i <- i + 1

  }

  tot_df <- do.call(rbind, tot_df) %>% select(-method)
  return(tot_df)

}

ModifyFisherOutput      <- function(df){

  ## Modify output from CalculateFisherEpitope-function

  # @ params:
  # df = df, from CalculateFisherEpitope -function

  df %>%
    mutate(p.adj = p.adjust(p.value, method = "BH"),
           sigf = ifelse(p.adj < 0.05, "sigf", "unsigf"),
           cluster_pred = cluster_predicted / c(cluster_predicted + cluster_unpredicted)) %>%
    arrange(p.adj)

}


RunFisherEpitope        <- function(viz_df, alternative = alternative){

  # @ params
  # input:
  # viz_df = df, that contains the TCRGP predictions and variables to count enrichment on
  # alternative = char, alternative for Fisher's exact test, 'less', 'greater' or 'two-sided'

  # output:
  # df, a tidy data-frame for erichment for pre-defined covariates


  ## The factors to consider when calculating Fisher's
  patients <- unique(viz_df$Patient)
  tissues  <- unique(viz_df$Sample.Type)
  stages   <- unique(viz_df$Stage)
  epitopes <- c("pol282_b", "pol387_b", "core169_b", "core195_b")
  # epitopes <- c("pol282_b", "pol387_b", "core169_b", "core195_b", "IPSINVHHY_trb", "TPRVTGGGAM_trb", "NLVPMVATV_trb", "GLCTLVAML_trb", "RAKFKQLL_trb", "YVLDHLIVV_trb", "GILGFVFTL_trb", "PKYVKQNTLKLAT_trb")

  vec2list <- function(vector){

    ## Helper function to translate vector into list, and the last element will be pooled from all the other elements

    vector_list <- list()
    for(i in 1:length(vector)){
      vector_list[[i]] <- as.character(vector[i])
    }

    vector_list[[c(length(vector)+1)]] <- as.character(vector_list)
    return(vector_list)
  }

  patient_list <- vec2list(vector = patients)
  tissue_list  <- vec2list(vector = tissues)
  stage_list   <- vec2list(vector = stages)
  epitope_list <- vec2list(vector = epitopes)

  ## Add I/II and IVB stage patients
  patient_list[[7]] <- c("P0205", "P0508", "P0407")
  patient_list[[8]] <- c("P1116", "P0322")

  epitope_list[[5]] <- c("hbv_reactive")

  fisher_list <- list()
  i <- 1

  for(patient_temp in patient_list){

    for(tissue_temp in tissue_list){

      for(stage_temp in stage_list){

        for(epitope_temp in epitope_list){

          fisher_list[[i]] <- CalculateFisherEpitope(viz_df = viz_df,
                                                     epitope_temp = epitope_temp,
                                                     patient_temp = patient_temp,
                                                     tissue_temp  = tissue_temp,
                                                     stage_temp   = stage_temp,
                                                     alternative  = alternative)
           i <- i + 1

        }
      }
    }
  }

  fisher_df <- do.call(rbind, fisher_list)
  fisher_df <- fisher_df[!duplicated(fisher_df), ]

  return(fisher_df)

}


getGini <- function(viz_df, covariate){

  ## Count clonality for each cluster

  # @ params:
  # viz_df    = df, that contains the TCRGP predictions and variables to count gini on
  # covariate = char, covariate to subset on

  ## Subset only to cluster

  cluster_gini <- NULL; i <- 1
  for(cluster in unique(viz_df[,covariate])){

    cluster_df             <- viz_df[viz_df[,covariate] == cluster,] %>% mutate(Clonotype.ID.x = as.factor(Clonotype.ID.x))
    cluster_df             <- cluster_df[!is.na(cluster_df$Clonotype.ID.x), ]
    cluster_gini[[i]]      <- data.frame(gini   = ineq::Gini(cluster_df$Clonotype.ID.x, na.rm = T),
                                         nTCR   = length(levels(cluster_df$Clonotype.ID.x)),
                                         nCells = nrow(cluster_df),
                                         factor = cluster)
    names(cluster_gini)[i] <- cluster

    i <- i + 1
  }

  gini_df <- do.call(rbind, cluster_gini) %>%
    as.data.frame() %>%
    mutate(test = covariate) %>% add_rownames(var = "covariate") %>% arrange(desc(gini))

  return(gini_df)

}



plotUMAP <- function(viz_df, cluster){

  viz_df_temp <- data.frame(viz_df, "seurat_cluster" = cluster)
  nClusters   <- unique(cluster) %>% length

  ## Visualise
  umap_mean <- data.frame(aggregate(UMAP1 ~ seurat_cluster, viz_df_temp, median), UMAP2 = aggregate(UMAP2 ~ seurat_cluster, viz_df_temp, median)[,2])


  ## Plot UMAPs with TCRGP predictions highlighted
  ggplot() +
    geom_point(data = viz_df_temp, aes(x = UMAP1, y = UMAP2, color = seurat_cluster), size = 0.8) +
    geom_point(data = subset(viz_df_temp, hbv_reactive == T), aes(x = UMAP1, y = UMAP2, fill = seurat_cluster), size = 2, shape = 21) +

    stat_ellipse(data = viz_df_temp, geom = "polygon", aes(x = UMAP1, y = UMAP2, color = seurat_cluster, fill = seurat_cluster), alpha = 0.1, lty = "dotted") +
    ggrepel::geom_label_repel(data = umap_mean,
                              aes(x = UMAP1, y = UMAP2, color = seurat_cluster, label = seurat_cluster),
                              size = 5, color = "black",
                              nudge_y       = 3,
                              segment.size  = 0.2,
                              segment.color = "grey50",
                              direction     = "x") +
    theme_void() +
    theme(legend.position = "none") +
    scale_color_manual(values = getPalette(nClusters)) +
    scale_fill_manual(values = getPalette(nClusters)) + labs()


}


plotHeatmap <- function(cells_to_use, show_heatmap_legend = F){

  ## Specific function for this project to plot heatmap with different cells

  # @ param:
  # cells_to_use: numeric vector containing the cells to use to plot heatmap

  require(ComplexHeatmap)
  require(circlize)

  naive_markers         <- c("TCF7", "SELL", "LEF1", "CCR7")
  cytotoxic_markers     <- c("IL2", "IFNG", "GZMA", "GZMB", "GZMM", "GZMH", "PRF1", "GNLY")
  costimulatory_markers <- c("CD28", "TNFRSF14", "ICOS", "TNFRSF9")
  inhibitory_markers    <- c("PDCD1", "TIGIT", "CTLA4", "LAG3", "HAVCR2")
  cd8_em                <- c("GZMK", "CXCR4", "CXCR3", "CD44")

  marker_list           <- list(naive_markers, cytotoxic_markers, costimulatory_markers, inhibitory_markers, cd8_em)
  marker_genes          <- do.call("c", marker_list)
  marker_amount         <- do.call("c", lapply(marker_list, length))
  split_known           <- rep(c("6 Naive", "5 Cytotoxic", "3 Costimulatory", "4 Inhibitory", "2 Effector \n memory"), marker_amount)

  # Count mean by cluster
  hcc_df                <- hcc_cd8@assays$RNA@data[marker_genes,cells_to_use] %>% as.matrix()
  hcc_df                <- rbind("Cluster" = as.factor(droplevels(Idents(hcc_cd8)[cells_to_use])), hcc_df)
  colnames(hcc_df)      <- NULL

  nLevels               <- as.factor(Idents(hcc_cd8)[cells_to_use]) %>% droplevels %>% levels() %>% length() 
  cluster_means_known   <- lapply(1:nLevels, function(df, i) rowMeans(df[ ,which(df[1,] == i)]),  df = hcc_df)
  cluster_means_known   <- do.call(rbind, cluster_means_known)

  ## Z-score
  cluster_means_known   <- cluster_means_known[,-1]
  cluster_means_known_z <- t(scale(cluster_means_known))

  ## Annotation
  clusters              <- as.factor(Idents(hcc_cd8)[cells_to_use]) %>% droplevels %>% levels()
  ha                    <- HeatmapAnnotation(names = anno_text(clusters, rot = 90))

  Heatmap(cluster_means_known_z,
          bottom_annotation = ha,
          col = colorRamp2(c(-2, 0, 4), c("dodgerblue", "white", "salmon")),
          cluster_columns = F,
          cluster_rows = F,
          row_names_gp = gpar(fontface = 3),
          row_names_side = c("right"),
          show_heatmap_legend = show_heatmap_legend,
          split = split_known,
          gap = unit(5, "mm"),
          bottom_annotation_height = unit(2.5, "cm"),
          name = "Z-score \n(of scaled log-counts)")

}
