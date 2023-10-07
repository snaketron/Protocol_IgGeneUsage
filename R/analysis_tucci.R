.libPaths(new = "/mnt/nfs/simo/rpack/")
require(IgGeneUsage)
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


ighj <- get(load("/mnt/nfs/simo/IgGeneUsage_manuscript/data_Tucci/ighj.RData"))
ighj <- ighj[ighj$productive==T & ighj$condition == "HD" & ighj$replicate == "1",]
ighj$condition <- ighj$cell
ighj <- ighj[, c("patient", "condition", "cell_count", "j_top_call_combi")]
colnames(ighj) <- c("sample_id", "condition", "gene_usage_count", "gene_name")
ighj$gap <- stringr::str_count(string = ighj$gene_name, pattern = '\\|')
ighj <- ighj[ighj$gap<=1,]
ighj$gap <- NULL


o <- IgGeneUsage::GU(ud = ighj,
                     mcmc_warmup = 2500,
                     mcmc_steps = 5000,
                     mcmc_chains = 5,
                     mcmc_cores = 5,
                     adapt_delta = 0.999,
                     max_treedepth = 13)

save(o, file = "/mnt/nfs/simo/IgGeneUsage_manuscript/data_Tucci/dgu_ighj_hd_rep1.RData")



.libPaths(new = "/mnt/nfs/simo/rpack/")
require(IgGeneUsage)
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


ighv <- get(load("/mnt/nfs/simo/IgGeneUsage_manuscript/data_Tucci/ighv.RData"))
ighv <- ighv[ighv$productive==T & ighv$condition == "HD" & ighv$replicate == "1",]
ighv$condition <- ighv$cell
ighv <- ighv[, c("replicate", "patient", "condition", "cell_count", "v_top_call_combi")]
colnames(ighv) <- c("replicate", "sample_id", "condition", "gene_usage_count", "gene_name")
ighv <- get_similar_genes(x = ighv)
ighv$gene_name <- ighv$gene_key
ighv <- get_similar_genes(x = ighv)

ighv <- aggregate(gene_usage_count~replicate+sample_id+condition+gene_key, 
                  data = ighv, FUN = sum)
ighv$gene_name <- ighv$gene_key
ighv$gene_key <- NULL

y <- ighv[ighv$gene_name %in% names(which(apply(
  X = table(ighv$gene_name, ighv$condition), MARGIN = 1, FUN = prod)==0)), ]
ighv <- ighv[!ighv$gene_name %in% y$gene_name, ]

# ighv <- ighv[ighv$replicate == 1, ]
# ighv <- ighv[ighv$replicate == 2 & ighv$sample_id != "HD6" , ]

ighv$replicate <- NULL
o <- IgGeneUsage::GU(ud = ighv,
                     mcmc_warmup = 2500,
                     mcmc_steps = 5000,
                     mcmc_chains = 5,
                     mcmc_cores = 5,
                     adapt_delta = 0.999,
                     max_treedepth = 13)

save(o, file = "/mnt/nfs/simo/IgGeneUsage_manuscript/data_Tucci/dgu_ighv_hd_rep1.RData")



ggplot(data = o$dgu_summary)+
  facet_wrap(facets = ~gene_name, scales = "free")+
  geom_hline(yintercept = 0)+
  geom_point(aes(x = gene_name, y = es_mean, col = c_name),
             position = position_dodge(width = 0.8))+
  geom_errorbar(aes(x = gene_name, y = es_mean, ymin = es_L, ymax = es_H, col = c_name),
                position = position_dodge(width = 0.8), width = 0.1)
