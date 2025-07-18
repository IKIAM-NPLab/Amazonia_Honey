
# Package for label PCA loading plot
library(ggrepel)

ggplot(scores, aes(PC1, PC2)) +
  geom_point(aes(shape = Species, color = Species), size = 3) +
  geom_segment(data = ei_compouds, aes(x = 0, y = 0, xend = PC1*100, yend = PC2*100),
               arrow = arrow(length = unit(0.3, "cm"), type = "open", angle = 25),
               size = 0.01, color = "darkblue") +
  geom_label_repel(data = ei_compouds, aes(label = meta_table$Metabolite,
                                           x = PC1*100, y = PC2*100),
                   box.padding = 0.37, label.padding = 0.22, label.r = 0.30, size = 2.5,
                   force = 4, max.overlaps = 50, min.segment.length = 0) +
  ggforce::geom_mark_ellipse(aes(filter = Species == "T. angustula"),
                             show.legend = FALSE, expand = unit(3, 'mm'),
                             tol = 0.01) +
  ggforce::geom_mark_ellipse(aes(filter = Species == "M. fasciculata"),
                             show.legend = FALSE, expand = unit(3.5, 'mm'),
                             tol = 0.01) +
  #ggforce::geom_mark_hull(aes(filter = Species == "M. fuscopilosa"),
  #                           show.legend = FALSE, expand = unit(3.5, 'mm'),
  #                           tol = 0.01) +
  geom_text(aes(label = Site), hjust = 0, nudge_x = 0.50,
            check_overlap = TRUE, show.legend = FALSE) +
  guides(x=guide_axis(title = "PC1 (20.14 %)"),
         y=guide_axis(title = "PC2 (17.62 %)")) +
  labs(shape = 'Bees species', color= 'Bees species') +
  theme_classic() +
  theme(legend.text = element_text(face="italic")) +
  theme(legend.position = c(0.120, 0.230),
        legend.background = element_rect(fill = "white", color = "black")) +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  geom_vline(xintercept = 0, linetype = "longdash", colour="gray") +
  geom_hline(yintercept = 0, linetype = "longdash", colour="gray")

