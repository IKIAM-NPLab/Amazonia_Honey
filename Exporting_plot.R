
# Exporting plot

## PCA of honey profile and antimicrobial activity

### Library loading
library(gridExtra)
### Figure matrix
figure_1 <- arrangeGrob(pca_plot,
                        corr_pca_plot,
                        load_pca,
                        corr_load_plot,
                        layout_matrix = rbind(c(1,   2),
                                              c(1,  2),
                                              rep(3, 2),
                                              rep(3, 2),
                                              rep(3, 2),
                                              rep(4, 2),
                                              rep(4, 2),
                                              rep(4, 2)))
### Adding label to the figures
figure_one <- ggpubr::as_ggplot(figure_1) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.5, 0, 0),
                  y = c(.99, .99, .75, .375))
### Exporting (*.pdf) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure1.pdf", plot = figure_one,
      width = 140, height = 180, units = "mm", dpi = 300, scale = 2.5)
### Exporting (*.png) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure1.png", plot = figure_one,
      width = 140, height = 180, units = "mm", dpi = 300, scale = 2.5)
### Exporting (*.jpg) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure1.jpg", plot = figure_one,
      width = 140, height = 180, units = "mm", dpi = 300, scale = 2.5)

## Heatmap with HCA and heatmap of the Pearson correlation

set.seed(1540)

# Add top anotation to HeatMap
top_info_ann_f2 <- HeatmapAnnotation(`Species` = hm_pdata$Species,
                                     col = list(`Species` = cols_species),
                                     annotation_name_side = "left",
                                     show_annotation_name = T,
                                     show_legend = F,
                                     border = T)

# Metabolomics heatmap
hm_plot_f2 <- Heatmap(hm_scl,
                      col = mycol,
                      border_gp = grid::gpar(col = "black", lty = 0.02),
                      rect_gp = grid::gpar(col = "black", lwd = 0.75),
                      clustering_distance_columns = "euclidean",
                      clustering_method_columns = "complete",
                      top_annotation = top_info_ann_f2,
                      column_names_gp = gpar(fontface = "italic"),
                      row_names_max_width = unit(10, "cm"),
                      #right_annotation = hm_row_ann,
                      show_heatmap_legend = F,
                      row_km = 3, column_km = 2,
                      row_title = c("a", "b", "c"))
# Metabolomics and correlation heatmap
hm_met_corr <- hm_plot_f2 + hm_corr
#Adding legends to heatmap
# Color scale legend
lgd1 <- Legend(col_fun = mycol,
               title = "Autoscaled abundance",
               direction = "horizontal" )
# Color scale legend
lgd1a <- Legend(col_fun = col_fun,
                title = "Correlation coefficient",
                direction = "horizontal" )
# Bees species legend
lgd2 <- Legend(labels = gt_render(c("*T. angustula*",
                                    "*M. fasciculata*",
                                    "*M. fuscopilosa*")),
               legend_gp = gpar(fill = cols_species),
               title = "Bees species", ncol = 1)
# Metabolite class Legend
lgd3 <- Legend(labels = c(unique(hm_fdata$classyfireR_Superclass)) ,
               legend_gp = gpar(fill = cols_metclass), 
               title = "Metabolite superclass", ncol = 2)
# Converting to ggplot
gg_f2 <- grid.grabExpr(draw(hm_met_corr))
gg_f2 <- ggpubr::as_ggplot(gg_f2)
# Legends
all_legends_f2 <- packLegend(lgd1, lgd2, lgd3, lgd1a, direction = "horizontal")
gg_legend_f2 <- grid.grabExpr(draw(all_legends_f2))
gg_legend_fn_f2 <- ggpubr::as_ggplot(gg_legend_f2)
# Heatmap plot
figure2 <- plot_grid(gg_legend_fn_f2,
                     gg_f2, ncol = 1,
                     rel_heights = c(0.055, 0.880))
### Exporting (*.pdf) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure2.pdf", plot = figure2,
       width = 160, height = 120, units = "mm", dpi = 300, scale = 2.5)
### Exporting (*.png) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure2.png", plot = figure2,
       width = 160, height = 120, units = "mm", dpi = 300, scale = 2.5)
### Exporting (*.jpg) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure2.jpg", plot = figure2,
       width = 160, height = 120, units = "mm", dpi = 300, scale = 2.5)

