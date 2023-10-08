.libPaths(new = "/mnt/nfs/simo/rpack/")
require(ggplot2)
require(ggforce)
require(patchwork)
require(rstan)

load("/mnt/nfs/simo/IgGeneUsage_manuscript/data_Tucci/fdata_5x.RData")

ighv_1_69 <- fdata_5x[fdata_5x$replicate == "1" & 
                        fdata_5x$productive == TRUE & 
                        fdata_5x$v_top_call_combi == "IGHV1-69x|" &
                        # fdata_5x$v_top_call_combi == "IGHV4-59x|" &
                        # fdata_5x$v_top_call_combi == "IGHV4-34x|" &
                        fdata_5x$patient %in% c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", "HD7"), ]

b <- fdata_5x[fdata_5x$replicate == "1" & 
                        fdata_5x$productive == TRUE & 
                        fdata_5x$v_top_call_combi != "IGHV1-69x|" &
                        fdata_5x$patient %in% c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", "HD7"), ]

f <- fdata_5x[fdata_5x$replicate == "1" & 
                fdata_5x$productive == TRUE & 
                fdata_5x$patient %in% c("HD1", "HD2", "HD3", "HD4", "HD5", "HD6", "HD7"), ]

rm(fdata_5x)


table(ighv_1_69$cell)
table(ighv_1_69$patient)
table(ighv_1_69$patient, ighv_1_69$cell)
table(ighv_1_69$patient, ighv_1_69$cell)/apply(table(ighv_1_69$patient, ighv_1_69$cell), 
                                               MARGIN = 1, FUN = sum)*100


# CDR3 length is not a marker
ggplot()+
  facet_wrap(facets = ~cell, nrow = 1)+
  geom_sina(data = b, aes(x = patient, y = junction_aa_length), col = "red", size = 0.5)+
  geom_sina(data = ighv_1_69, aes(x = patient, y = junction_aa_length), size = 1)
  
ggplot()+
  facet_wrap(facets = ~cell, nrow = 1)+
  geom_point(data = aggregate(junction_aa_length~patient+cell, data = ighv_1_69, FUN = mean),
             aes(x = patient, y = junction_aa_length), col = "black")+
  geom_point(data = aggregate(junction_aa_length~patient+cell, data = b, FUN = mean),
             aes(x = patient, y = junction_aa_length), col = "red")


ggplot()+
  facet_wrap(facets = ~cell, nrow = 1)+
  geom_sina(data = aggregate(junction_aa_length~patient+cell+v_top_call_combi, data = b, FUN = mean),
             aes(x = patient, y = junction_aa_length), col = "red")+
  geom_point(data = aggregate(junction_aa_length~patient+cell, data = ighv_1_69, FUN = mean),
             aes(x = patient, y = junction_aa_length), col = "black", size = 3)

aggregate(junction_aa_length~patient+cell, data = ighv_1_69, FUN = mean)


# ighv_1_69$charge <- alakazam::charge(ighv_1_69$junction_aa, pH = 7.4, pK = NULL, normalize = TRUE)
ighv_1_69$charge <- Peptides::charge(seq = ighv_1_69$junction_aa)

# CDR3 charge is not a marker
ggplot(data = ighv_1_69)+
  facet_wrap(facets = ~cell, nrow = 1)+
  geom_sina(aes(x = patient, y = charge), size = 1)+
  geom_point(data = aggregate(charge~patient+cell, data = ighv_1_69, FUN = mean),
             aes(x = patient, y = charge), size = 2, col = "red")

  

require(rstan)
nb <- rstan::stan_model(file = "../IgGeneUsage_manuscript/data_Tucci/nb.stan")

get_model_len <- function(x, d, m) {
  cat(x, "\n")
  d <- d[d$v_top_call_combi == x, ]
  d$c <- as.numeric(as.factor(d$cell))
  if(max(d$c) == 4) {
    f <- rstan::sampling(object = m,
                         data = list(N = nrow(d),
                                     y = d$junction_aa_length,
                                     x = as.numeric(as.factor(d$cell))),
                         chains = 4,
                         cores = 1,
                         iter = 2000,
                         warmup = 1000)
    
    d <- d[duplicated(d[, c("cell", "c", "v_top_call_combi")])==F,]
    d <- d[, c("cell", "c", "v_top_call_combi")]
    
    summary_mu <- data.frame(summary(f, par = "mu")$summary)
    summary_mu$c <- as.numeric(gsub(pattern = "mu\\[|\\]", replacement = '', x = rownames(summary_mu)))
    summary_mu <- summary_mu[, c("mean", "X2.5.", "X97.5.", "c")]
    colnames(summary_mu) <- paste0("mu_", colnames(summary_mu))
    
    summary_sigma <- data.frame(summary(f, par = "sigma")$summary)
    summary_sigma$c <- as.numeric(gsub(pattern = "sigma\\[|\\]", replacement = '', x = rownames(summary_sigma)))
    summary_sigma <- summary_sigma[, c("mean", "X2.5.", "X97.5.", "c")]
    colnames(summary_sigma) <- paste0("sigma_", colnames(summary_sigma))
    
    summary_mu <- merge(x = summary_mu, y = summary_sigma, by.x = "mu_c", by.y = "sigma_c")
    summary_mu <- merge(x = summary_mu, y = d, by.x = "mu_c", by.y = "c")
    
    return(summary_mu)
  }
  
  return(NULL)
}

out <- parallel::mclapply(X = unique(f$v_top_call_combi), 
                          FUN = get_model_len, 
                          m = nb, 
                          d = f,
                          mc.cores = 20)
save(out, file = "../IgGeneUsage_manuscript/data_Tucci/nb_lengths.RData")
o <- do.call(rbind, out)


ggplot(data = o)+
  facet_wrap(facets = ~v_top_call_combi, scales = "free_y")+
  geom_point(aes(x = cell, y = mu_mean))+
  geom_errorbar(aes(x = cell, y = mu_mean, ymin = mu_X2.5., ymax = mu_X97.5.),width = 0.1)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
