
## Visualise
message("Visualize results...")




## Plot UMAPs with TCRGP predictions highlighted
p <- plotUMAP(viz_df, viz_df$cluster)
ggsave(plot = p, paste0(output_dir, "umap_tcrgp_prediction.pdf"), width = 8, height = 6)




## Plot barplot
viz_df$meta_epitope <- ifelse(viz_df$epitope_spec == "None", viz_df$viral_epitope_spec, viz_df$epitope_spec)
bar_df              <- viz_df %>% group_by(cluster, meta_epitope) %>% summarise(n = n()) %>% mutate(prop = n/sum(n))

bar_df$meta_epitope <- as.factor(bar_df$meta_epitope)
bar_df$meta_epitope <- relevel(bar_df$meta_epitope, ref = "None")

levels(bar_df$meta_epitope) <- c("No epitope",
                                 expression("HBV"["core"~169]), expression("HBV"["core"~195]),
                                 expression("CMV pp65"["123"]), expression("CMV pp65"["495"]), expression("Inf HA"["306"]),
                                 expression("HBV"["pol"~282]), expression("HBV"["pol"~387]),
                                 expression("EBV BZLF1"["190"]),  expression("CMV pp65"["417"]),  expression("EBV BRLF1"["109"]))

bar_df$cluster      <- factor(as.character(bar_df$cluster), levels = levels(bar_df$cluster)[order(levels(bar_df$cluster))])
bar_df$meta_epitope <- factor(as.character(bar_df$meta_epitope), levels = levels(bar_df$meta_epitope)[order(levels(bar_df$meta_epitope))])
bar_labels <- levels(bar_df$meta_epitope)

bar_labels <- c(expression("CMV pp65"["123"]), expression("CMV pp65"["417"]), expression("CMV pp65"["495"]),
                expression("EBV BRLF1"["109"]), expression("EBV BZLF1"["190"]),
                expression("HBV"["core"~169]), expression("HBV"["core"~195]), expression("HBV"["pol"~282]), expression("HBV"["pol"~387]), expression("IAV HA"["306"]), "No epitope")


# bar_cols <- c("lightgrey", brewer.pal(4, "OrRd"), brewer.pal(6, "Blues"))
bar_cols <- c(brewer.pal(3, "Greens"), brewer.pal(2, "Blues")[1:2], brewer.pal(4, "OrRd"), "rosybrown2", "lightgrey")

p <- ggplot(bar_df, aes(cluster, prop, fill = meta_epitope, label=n)) +
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +

  scale_fill_manual(values = bar_cols, labels = bar_labels) +
  labs(x = "", y = "fraction of cells", fill = "TCRGP prediction") +
  theme_classic() + theme(legend.text.align = 0) +
  theme(text = element_text(size = 12),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12))
ggsave(plot = p, paste0(output_dir, "bar_tcrgp_prediction.pdf"), width = 8, height = 4)




## Plot heatmaps
exhausted_cells <- viz_df$cluster %in% c("exhausted 1", "exhausted 2", "exhausted 3")

pdf(paste0(output_dir, "heatmap_epitope_scaled.pdf"), width = 3, height = 10)
plotHeatmap(cells_to_use = exhausted_cells) %>% print
dev.off()

pdf(paste0(output_dir, "heatmap_epitope_for_scale.pdf"), width = 3, height = 10)
plotHeatmap(cells_to_use = exhausted_cells, show_heatmap_legend = T) %>% print
dev.off()
