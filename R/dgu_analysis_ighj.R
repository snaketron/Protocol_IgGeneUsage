.libPaths(new = "/mnt/nfs/simo/rpack/")
require(IgGeneUsage)
require(knitr)
require(ggplot2)
require(ggforce)
require(patchwork)


load("../IgGeneUsage_manuscript/data_Tucci/ighj.RData")
d <- ighj[ighj$replicate == "1" & ighj$productive == TRUE, ]
d <- d[,c("condition", "cell", "cell_count", "patient", "j_top_call_combi")]
d$gene_name <- d$j_top_call_combi
d$j_top_call_combi <- NULL
d <- d[which(stringr::str_count(string = d$gene_name, pattern = "\\|") == 1), ]


# keep only HCV
hcv <- d[d$condition != "HD", ]
hcv$condition <- NULL
colnames(hcv) <- c("condition", "gene_usage_count", "sample_id", "gene_name")

hcv <- aggregate(gene_usage_count~sample_id+condition+gene_name, data = hcv, FUN = sum)


# keep only HD
hd <- d[d$condition == "HD", ]
hd$condition <- NULL
colnames(hd) <- c("condition", "gene_usage_count", "sample_id", "gene_name")

hd <- aggregate(gene_usage_count~sample_id+condition+gene_name, data = hd, FUN = sum)


# keep HD and HCV
hd_hcv <- d
colnames(hd_hcv) <- c("condition", "cell", "gene_usage_count", "sample_id", "gene_name")

hd_hcv <- aggregate(gene_usage_count~sample_id+condition+gene_name+cell, data = hd_hcv, FUN = sum)






# data_ighj <- list(hd = hd, hcv = hcv, hd_hcv = hd_hcv)
# save(data_ighj, file = "dgu_Tucci/data_ighj.RData")
data_ighj <- get(load("dgu_Tucci/data_ighj.RData"))


dgu_hd <- IgGeneUsage::DGU(ud = data_ighj$hd,
                           mcmc_warmup = 1500,
                           mcmc_steps = 3000,
                           mcmc_chains = 5,
                           mcmc_cores = 5,
                           adapt_delta = 0.99,
                           max_treedepth = 13)
save(dgu_hd, file = "dgu_Tucci/DGU_hd_ighj_rep1.RData")




dgu_hcv <- IgGeneUsage::DGU(ud = data_ighj$hcv,
                           mcmc_warmup = 1500,
                           mcmc_steps = 3000,
                           mcmc_chains = 5,
                           mcmc_cores = 5,
                           adapt_delta = 0.99,
                           max_treedepth = 13)
save(dgu_hcv, file = "dgu_Tucci/DGU_hcv_ighj_rep1.RData")


