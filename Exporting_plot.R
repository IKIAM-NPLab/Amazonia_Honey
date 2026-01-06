
# Exporting plot

## PCA of honey profile and antimicrobial activity

### Library loading
library(gridExtra)
library(cowplot)
### Figure matrix
figure_1 <- arrangeGrob(pca_plot,
                        corr_pca_plot,
                        #load_pca,
                        corr_load_plot,
                        layout_matrix = rbind(c(1,   2),
                                              c(1,  2),
                                              rep(3, 2),
                                              rep(3, 2),
                                              rep(3, 2)#,
                                              #rep(4, 2),
                                              #rep(4, 2),
                                              #rep(4, 2)
                                              ))
### Adding label to the figures
figure_one <- ggpubr::as_ggplot(figure_1) +
  draw_plot_label(label = LETTERS[1:3],
                  x = c(0, 0.5, 0),
                  y = c(.99, .99, .60))
### Exporting (*.pdf) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure1.pdf", plot = figure_one,
      width = 140, height = 100, units = "mm", dpi = 300, scale = 2.3)
### Exporting (*.png) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure1.png", plot = figure_one,
      width = 140, height = 100, units = "mm", dpi = 300, scale = 2.3)
### Exporting (*.jpg) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure1.jpg", plot = figure_one,
      width = 140, height = 100, units = "mm", dpi = 300, scale = 2.3)

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
                      row_names_max_width = unit(8, "cm"),
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
                                    "*M. grandis*",
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
                     rel_heights = c(1, 9))
### Exporting (*.pdf) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure2.pdf", plot = figure2,
       width = 140, height = 120, units = "mm", dpi = 300, scale = 2.0)
### Exporting (*.png) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure2.png", plot = figure2,
       width = 140, height = 120, units = "mm", dpi = 300, scale = 2.0)
### Exporting (*.jpg) file
ggsave(filename = "Result/notame_Result/HS_GCMS/Figure2.jpg", plot = figure2,
       width = 140, height = 120, units = "mm", dpi = 300, scale = 2.0)

## Antimicrobial activity

# Avoid x labels
dn_ec100_nxl <- dn_ec100 + theme(axis.text.x = element_blank())
dn_sa100_nxl <- dn_sa100 + theme(axis.text.x = element_blank())

### Figure matrix
figure_a3 <- arrangeGrob(dn_ec100_nxl,
                         dn_sa100_nxl,
                         dn_kp100,
                         tk_pm100,
                         layout_matrix = rbind(c(1, 2),
                                               c(3, 4)),
                         heights=c(2.9, 4))
### Adding label to the figures
figure_athree <- ggpubr::as_ggplot(figure_a3) +
  draw_plot_label(label = LETTERS[1:4],
                  x = c(0, 0.5, 0, 0.5),
                  y = c(.99, .99, .57, .57))
### Exporting (*.pdf) file
ggsave(filename = "Result/Antimicrobial/Figure_A3.pdf", plot = figure_athree,
       width = 110, height = 80, units = "mm", dpi = 300, scale = 2.5)
### Exporting (*.png) file
ggsave(filename = "Result/Antimicrobial/Figure_A3.png", plot = figure_athree,
       width = 110, height = 80, units = "mm", dpi = 300, scale = 2.5)
### Exporting (*.jpg) file
ggsave(filename = "Result/Antimicrobial/Figure_A3.jpg", plot = figure_athree,
       width = 110, height = 80, units = "mm", dpi = 300, scale = 2.5)

## PCA biplot

# Package for label the loading of the PCA biplot
library(ggrepel)
# PCA biplot result
fig_a4 <- ggplot(scores, aes(PC1, PC2)) +
  geom_point(aes(shape = Species, color = Species), size = 3) +
  scale_color_manual(values=c("#7CAE00",
                              "#F8766D",
                              "#00BFC4",
                              "#E76BF3")) +
  scale_shape_manual(values=c(17, 16, 15, 3)) +
  geom_segment(data = ei_compouds, aes(x = 0, y = 0, xend = PC1*100, yend = PC2*100),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed", angle = 25),
               size = 0.01, color = "darkblue") +
  geom_label_repel(data = ei_compouds, aes(label = meta_table$Metabolite,
                                           x = PC1*100, y = PC2*100),
                   box.padding = 0.37, label.padding = 0.22, label.r = 0.30, size = 3,
                   force = 4, max.overlaps = 50, min.segment.length = 0) +
  ggforce::geom_mark_ellipse(aes(filter = Species == "T. angustula", color = Species),
                             show.legend = FALSE, expand = unit(3, 'mm'),
                             tol = 0.01) +
  ggforce::geom_mark_ellipse(aes(filter = Species == "M. grandis", color = Species),
                             show.legend = FALSE, expand = unit(8.2, 'mm'),
                             tol = 0.01) +
  guides(x=guide_axis(title = "Loadings on PC1 (20.14 %)"),
         y=guide_axis(title = "Loadings on PC2 (17.62 %)")) +
  labs(shape = 'Bees species', color= 'Bees species') +
  theme_classic() +
  theme(legend.text = element_text(face="italic")) +
  theme(legend.position = c(0.060, 0.130),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")
### Exporting (*.pdf) file
ggsave(filename = "Result/notame_Result/HS_GCMS/FigureA4.pdf", plot = fig_a4,
       width = 168, height = 84, units = "mm", dpi = 300, scale = 2.0)
### Exporting (*.png) file
ggsave(filename = "Result/notame_Result/HS_GCMS/FigureA4.png", plot = fig_a4,
       width = 168, height = 84, units = "mm", dpi = 300, scale = 2.0)
### Exporting (*.jpg) file
ggsave(filename = "Result/notame_Result/HS_GCMS/FigureA4.jpg", plot = fig_a4,
       width = 168, height = 84, units = "mm", dpi = 300, scale = 2.0)
