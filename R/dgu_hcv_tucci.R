# .libPaths(new = "/mnt/nfs/simo/rpack/")
require(IgGeneUsage)
require(knitr)
require(ggplot2)
require(ggforce)
require(patchwork)


get_similar_genes <- function(x) {
  
  get_genes <- function(x) {
    return(unlist(strsplit(x = x, split = '\\|')))
  }
  
  get_matching_genes <- function(x, k) {
    # browser()
    ks <- unlist(strsplit(x = k, split = '\\|'))
    
    i <- any(x %in% ks)
    if(i) {
      return(x)
    }
    return(NULL)
  }
  
  gs <- lapply(X = x$gene_name, FUN = get_genes)
  
  x$gene_key <- NA
  for(i in 1:nrow(x)) {
    h <- lapply(X = gs, FUN = get_matching_genes, k = x$gene_name[i])
    x$gene_key[i] <- paste0(sort(unique(unlist(h))), collapse = '|')
  }
  return(x)
}

load("../IgGeneUsage_manuscript/data_Tucci/ighv.RData")
hd <- ighv[ighv$replicate == "1" & ighv$productive == TRUE, ]
hd <- hd[,c("condition", "cell", "cell_count", "patient", "v_top_call_combi")]
hd$gene_name <- hd$v_top_call_combi
hd$v_top_call_combi <- NULL

# get groups of genes
hd <- get_similar_genes(x = hd)
hd$gene_name <- hd$gene_key
hd <- get_similar_genes(x = hd)

# keep only HCV
hd <- hd[hd$condition != "HD", ]
hd$condition <- NULL
hd$gene_name <- hd$gene_key
hd$gene_key <- NULL

colnames(hd) <- c("condition", "gene_usage_count", "sample_id", "gene_name")

hd <- aggregate(gene_usage_count~sample_id+condition+gene_name, data = hd, FUN = sum)

y <- hd[hd$gene_name %in% names(which(apply(
  X = table(hd$gene_name, hd$condition), MARGIN = 1, FUN = prod)==0)), ]
hd <- hd[!hd$gene_name %in% y$gene_name, ]


dgu <- IgGeneUsage::DGU(ud = hd,
                        mcmc_warmup = 3000,
                        mcmc_steps = 8000,
                        mcmc_chains = 5,
                        mcmc_cores = 5,
                        adapt_delta = 0.99,
                        max_treedepth = 12)
save(dgu, file = "/mnt/nfs/simo/IgGeneUsage_manuscript/data_Tucci/DGU_hcv_rep1.RData")




require(IgGeneUsage)
require(knitr)
require(ggplot2)
require(ggforce)
require(patchwork)
dgu <- get(load("~/Desktop/work/R/IgGeneUsage_manuscript/data_Tucci/DGU_hcv_rep1.RData"))

rstan::check_hmc_diagnostics(object = dgu$glm)
rstan::stan_rhat(object = dgu$glm)
rstan::stan_ess(object = dgu$glm)


x <- dgu$dgu_summary
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


g2 <- ggplot(data = x)+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = gene_name, y = mult*es_mean, col = contrast),
             position = position_dodge(width = 0.6), size = 1)+
  geom_errorbar(aes(x = gene_name, y = mult*es_mean, col = contrast, 
                    ymin = mult*es_L, ymax = mult*es_H), width = 0,
                position = position_dodge(width = 0.6))+
  theme_bw()+
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(label = expression(gamma))+
  xlab(label = "Ig gene")



g1 <- ggplot(data = dgu$gu_summary)+
  geom_point(aes(x = gene_name, y = prob_mean, col = condition),
             position = position_dodge(width = 0.7), size = 1)+
  geom_errorbar(aes(x = gene_name, y = prob_mean, col = condition, 
                    ymin = prob_L, ymax = prob_H), width = 0.4,
                position = position_dodge(width = 0.7))+
  theme_bw(base_size = 10)+
  theme(legend.position = "top", 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab(label = expression(alpha))+
  xlab(label = "Ig gene")

g <- (g1/g2)
g

ggsave(filename = "../IgGeneUsage_manuscript/data_Tucci/HCV_GU.pdf", 
       plot = g1,
       device = "pdf",
       width = 8,
       height = 4)




x$es_mean_adj <- x$es_mean*x$mult
pdf(file = "../IgGeneUsage_manuscript/data_Tucci/tree_HCV.pdf",
    width = 9, height = 5.5)
w <- reshape2::acast(data = x, formula = gene_name~contrast, value.var = "es_mean_adj")
w <- plot(hclust(method = "ward.D", dist(w)), cex = 0.9)
dev.off()



g3 <- ggplot(data = x[x$gene_name %in% c("3-72", "6-1", "3-7", "3-74"),])+
  facet_wrap(facets = ~gene_name, scales = "free_x", ncol = 5)+
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




g4 <- ggplot(data = x[x$gene_name %in% c("1-69", "4-34", "1-24", "2-26", "1-58", "1-18", "1-45"),])+
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

g <- (g3|g4)+patchwork::plot_layout(widths = c(4, 7))
g
ggsave(filename = "../IgGeneUsage_manuscript/data_Tucci/HCV_DGUs_2.pdf", 
       plot = g,
       device = "pdf",
       width = 8,
       height = 2.35)

