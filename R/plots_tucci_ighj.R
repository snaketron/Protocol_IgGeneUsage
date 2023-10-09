.libPaths(new = "/mnt/nfs/simo/rpack/")
dir.create("results_Tucci_ighj/")
require(IgGeneUsage)
require(knitr)
require(ggplot2)
require(ggforce)
require(patchwork)

dgu_hd <- get(load("dgu_Tucci/DGU_hd_ighj_rep1.RData"))


# HD
x <- dgu_hd$dgu_summary
y <- do.call(rbind, strsplit(x = x$contrast, split = "\\-vs\\-"))
x$c1 <- y[,1]
x$c2 <- y[,2]
x$g1 <- y[,1]
x$g2 <- y[,2]
x$g1[!x$g1 %in% c("CS", "NCS")] <- "N"
x$g1[x$g1 %in% c("CS", "NCS")] <- "M"
x$g2[!x$g2 %in% c("CS", "NCS")] <- "N"
x$g2[x$g2 %in% c("CS", "NCS")] <- "M"
rm(y)

x$mult <- 1
x[x$g1 == "N" & x$g2 == "M", "mult"]<--1
x$contrast <- ifelse(test = x$mult == 1, 
                     yes = paste0(x$c1, "-vs-", x$c2), 
                     no = paste0(x$c2, "-vs-", x$c1))
x_hd <- x
rm(x)



plot_dgu_hd <- ggplot(data = x_hd)+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = gene_name, y = mult*es_mean, col = contrast),
             position = position_dodge(width = 0.6), size = 1)+
  geom_errorbar(aes(x = gene_name, y = mult*es_mean, col = contrast, 
                    ymin = mult*es_L, ymax = mult*es_H), width = 0,
                position = position_dodge(width = 0.6))+
  theme_bw(base_size = 10)+
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(label = expression(gamma))+
  xlab(label = "Ig gene")+
  scale_color_discrete(name = '')+
  guides(colour = guide_legend(nrow = 1))




plot_gu_hd <- ggplot(data = dgu_hd$gu_summary)+
  geom_point(aes(x = gene_name, y = prob_mean, col = condition),
             position = position_dodge(width = 0.7), size = .8)+
  geom_errorbar(aes(x = gene_name, y = prob_mean, col = condition, 
                    ymin = prob_L, ymax = prob_H), width = 0.2,
                position = position_dodge(width = 0.7))+
  theme_bw(base_size = 10)+
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(label = expression(alpha))+
  xlab(label = "Ig gene")


ggsave(filename = "results_Tucci_ighj/gu_hd.pdf", 
       plot = plot_gu_hd,
       device = "pdf",
       width = 2.5,
       height = 2.75)






x_hd$es_mean_adj <- x_hd$es_mean*x_hd$mult
pdf(file = "results_Tucci_ighj/tree_hd.pdf",
    width = 4.4, height = 5)
w <- reshape2::acast(data = x_hd, 
                     formula = gene_name~contrast, 
                     value.var = "es_mean_adj")
w <- plot(hclust(method = "ward.D", dist(w)), cex = 0.7)
dev.off()



gn <- ggplot(data = x_hd[x_hd$gene_name %in% c("IGHJ6|"),])+
  facet_wrap(facets = ~gene_name, scales = "free_x", nrow = 1)+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = '', y = mult*es_mean, col = contrast),
             position = position_dodge(width = 1), size = 1)+
  geom_errorbar(aes(x = '', y = mult*es_mean, col = contrast, 
                    ymin = mult*es_L, ymax = mult*es_H), width = 0,
                position = position_dodge(width = 1))+
  theme_bw(base_size = 10)+
  theme(legend.position = "top")+
  ylab(label = expression(gamma))+
  xlab(label = "Ig gene")



gp <- ggplot(data = x_hd[x_hd$gene_name %in% c("IGHJ1|", "IGHJ2|", "IGHJ4|"),])+
  facet_wrap(facets = ~gene_name, scales = "free_x", nrow = 1)+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = '', y = mult*es_mean, col = contrast),
             position = position_dodge(width = 1), size = 1)+
  geom_errorbar(aes(x = '', y = mult*es_mean, col = contrast, 
                    ymin = mult*es_L, ymax = mult*es_H), width = 0,
                position = position_dodge(width = 1))+
  theme_bw(base_size = 10)+
  theme(legend.position = "top")+
  ylab(label = expression(gamma))+
  xlab(label = "Ig gene")

g <- (gn|gp)+patchwork::plot_layout(widths = c(1, 3))

ggsave(filename = "results_Tucci_ighj/tree_clusters.pdf", 
       plot = g,
       device = "pdf",
       width = 5,
       height = 2.35)


heatmap <- ggplot(data = x_hd)+
  geom_tile(aes(x = contrast, y = gene_name, fill = es_mean_adj), col = "white")+
  theme_bw(base_size = 10)+
  # theme(legend.position = "top")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab(label = "Contrasts")+
  ylab(label = "Ig gene")+
  scale_fill_distiller(name = expression(gamma), palette = "Spectral")+
  coord_flip()
  
ggsave(filename = "results_Tucci_ighj/tree_example.pdf", 
       plot = heatmap,
       device = "pdf",
       width = 3.2,
       height = 1.75)


