.libPaths(new = "/mnt/nfs/simo/rpack/")
require(IgGeneUsage)
require(knitr)
require(ggplot2)
require(ggforce)
require(patchwork)


group_genes <- function(x) {
  
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
d <- ighv[ighv$replicate == "1" & ighv$productive == TRUE, ]
d <- d[,c("condition", "cell", "cell_count", "patient", "v_top_call_combi")]
d$gene_name <- d$v_top_call_combi
d$v_top_call_combi <- NULL

# get groups of genes
d <- group_genes(x = d)
d$gene_name <- d$gene_key
d <- group_genes(x = d)

# keep only HCV
hcv <- d[d$condition != "HD", ]
hcv$condition <- NULL
hcv$gene_name <- hcv$gene_key
hcv$gene_key <- NULL
colnames(hcv) <- c("condition", "gene_usage_count", "sample_id", "gene_name")

hcv <- aggregate(gene_usage_count~sample_id+condition+gene_name, data = hcv, FUN = sum)

# y <- hcv[hcv$gene_name %in% names(which(apply(
#   X = table(hcv$gene_name, hcv$condition), MARGIN = 1, FUN = prod)==0)), ]
# hcv <- hcv[!hcv$gene_name %in% y$gene_name, ]


# keep only HD
hd <- d[d$condition == "HD", ]
hd$condition <- NULL
hd$gene_name <- hd$gene_key
hd$gene_key <- NULL
colnames(hd) <- c("condition", "gene_usage_count", "sample_id", "gene_name")

hd <- aggregate(gene_usage_count~sample_id+condition+gene_name, data = hd, FUN = sum)


# table(hd$gene_name,hd$sample_id)
# apply(X = table(hd$sample_id, hd$gene_name), MARGIN = 1, FUN = sum)
# y <- hd[hd$gene_name %in% names(which(apply(
#   X = table(hd$gene_name, hd$condition), MARGIN = 1, FUN = prod)==0)), ]
# hd <- hd[!hd$gene_name %in% y$gene_name, ]


dgu_hd <- IgGeneUsage::DGU(ud = hd,
                           mcmc_warmup = 3000,
                           mcmc_steps = 8000,
                           mcmc_chains = 5,
                           mcmc_cores = 5,
                           adapt_delta = 0.999,
                           max_treedepth = 13)
save(dgu_hd, file = "dgu_Tucci/DGU_hd_rep1.RData")




dgu_hcv <- IgGeneUsage::DGU(ud = hcv,
                           mcmc_warmup = 3000,
                           mcmc_steps = 8000,
                           mcmc_chains = 5,
                           mcmc_cores = 5,
                           adapt_delta = 0.999,
                           max_treedepth = 13)
save(dgu_hcv, file = "dgu_Tucci/DGU_hcv_rep1.RData")

cat("DONE")
