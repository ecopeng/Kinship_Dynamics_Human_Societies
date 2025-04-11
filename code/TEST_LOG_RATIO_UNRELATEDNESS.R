rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### load the predicted and observed unrelatedness
################################################################################
load('./data/relatedness/PREDICTIONS/log_ratio_UR_f2m_PRED.rdata')
load('./data/relatedness/OBSERVATIONS/log_ratio_UR_f2m_OBSE.rdata')

PRED <- log_ratio_UR_f2m_PRED
OBSE <- log_ratio_UR_f2m_OBSE
OBSE <- na.omit(OBSE)

MSD <- function(df_DATA, df_PRED) {
  df_DATA <- na.omit(df_DATA)
  idx_in_pred <- match(df_DATA$age, df_PRED$age)
  df_DATA <- df_DATA[!is.na(idx_in_pred), ]
  idx_in_pred <- na.omit(idx_in_pred)
  df_PRED <- df_PRED[idx_in_pred, ]
  DAT <- df_DATA$log_ur_f2m
  PRE <- df_PRED$log_ur_f2m
  return(sum((PRE - DAT)^2) / length(DAT))
}

vMSD <- function(df_dat, df_mod, society, per_soc_nmnf) {
  vmsd <- c()
  for(p in 1:length(per_soc_nmnf)) {
    df_O <- df_dat[which(df_dat$society == society[p]), ]
    soc_nm_nf <- strsplit(per_soc_nmnf[p], '\\+')[[1]]
    df_P <- df_mod[which(df_mod$society == soc_nm_nf[1] &
                         df_mod$Nm == as.numeric(soc_nm_nf[2]) &
                         df_mod$Nf == as.numeric(soc_nm_nf[3])), ]
    vmsd <- c(vmsd, MSD(df_O, df_P))
  }
  return(vmsd)
}

################################################################################
##### find the best fit Nm and Nf in the model (unrelatedness)
##### i.e., the best-fit scenarios are the same as the best-fit relatedness
################################################################################
society <- levels(OBSE$society)

load('./data/relatedness/PREDICTIONS/NmNf_relatedness_best_fit.rdata')
summary_t_msd_r_fm$soc_data <- factor(summary_t_msd_r_fm$soc_data, levels = society)
summary_t_msd_r_fm <- summary_t_msd_r_fm %>% arrange(soc_data)
if(!identical(society, as.character(summary_t_msd_r_fm$soc_data))) {stop('factor levels for *society* inconsistent!')}
if(!identical(society, levels(summary_t_msd_r_fm$soc_data))) {stop('factor levels for *society* inconsistent!')}

log_ratio_unrelatedness_best_fit <- data.frame()
for(u in 1:nrow(summary_t_msd_r_fm)) {
  log_ratio_unrelatedness_best_fit <- rbind(log_ratio_unrelatedness_best_fit,
               PRED[which(PRED$society == summary_t_msd_r_fm[u, ]$soc_pred &
                          PRED$Nm == summary_t_msd_r_fm[u, ]$Nm &
                          PRED$Nf == summary_t_msd_r_fm[u, ]$Nf), ])
}
save(log_ratio_unrelatedness_best_fit, file = './data/relatedness/PREDICTIONS/log_ratio_unrelatedness_best_fit.rdata')

################################################################################
##### test model predictions: permutation-based approach
################################################################################
v_soc_NmNf <- paste(summary_t_msd_r_fm$soc_pred,  # soc id
                    summary_t_msd_r_fm$Nm,        # Nm
                    summary_t_msd_r_fm$Nf,        # Nf
                    sep = '+')
if(!identical(substr(society, 4, 6), substr(v_soc_NmNf, 1, 3))) {stop('factor levels for *society* inconsistent!')}

tab_soc_NmNf <- table(v_soc_NmNf)
mat_per_soc_nmnf_soc <- t(unique(permut(v_soc_NmNf)))

deno <- 1
for(u in 1:length(tab_soc_NmNf)) {deno <- deno * factorial(as.numeric(tab_soc_NmNf[u]))}
len <- factorial(sum(tab_soc_NmNf)) / deno

MSDs <- matrix(nrow = length(society), ncol = ncol(mat_per_soc_nmnf_soc))
for(u in 1:ncol(mat_per_soc_nmnf_soc)) {
  MSDs[, u] <- vMSD(OBSE, log_ratio_unrelatedness_best_fit, society, mat_per_soc_nmnf_soc[, u])
}

distr_PER <- data.frame(sMSD = colSums(MSDs))
quant_obs_PER <- sum(vMSD(OBSE, log_ratio_unrelatedness_best_fit, society, v_soc_NmNf))
quant_cri_PER <- quantile(distr_PER, probs = .05, na.rm = T)
cat('PERMUTATION TEST:', '\n', 'actual quantile: ', quant_obs_PER, '\n', 'critical quantile: ', quant_cri_PER, '\n')

log_ratio_unrelatedness_permutation_test_results <- list(distr = distr_PER, observed_quantile = quant_obs_PER, critical_quantile = quant_cri_PER)
save(log_ratio_unrelatedness_permutation_test_results, file = './data/analyses/log_ratio_unrelatedness_permutation_test_results.rdata')

################################################################################
##### test model predictions: randomization-based approach
################################################################################
p <- unname(tab_soc_NmNf) / sum(unname(tab_soc_NmNf))
# the number of random samples = the number of unique permutations (len)
MSDs <- matrix(nrow = sum(unname(tab_soc_NmNf)), ncol = len)
for(u in 1:len) {
  soc_nmnf_rnd <- sample(names(tab_soc_NmNf),
                         sum(unname(tab_soc_NmNf)),
                         replace = T, prob = p)
  MSDs[, u] <- vMSD(OBSE, log_ratio_unrelatedness_best_fit, society, soc_nmnf_rnd)
}

distr_RND <- data.frame(sMSD = colSums(MSDs))
quant_obs_RND <- sum(vMSD(OBSE, log_ratio_unrelatedness_best_fit, society, v_soc_NmNf))
quant_cri_RND <- quantile(distr_RND, probs = .05, na.rm = T)
cat('RANDOMIZATION TEST:', '\n', 'actual quantile: ', quant_obs_RND, '\n', 'critical quantile: ', quant_cri_RND, '\n')

log_ratio_unrelatedness_randomization_test_results <- list(distr = distr_RND, observed_quantile = quant_obs_RND, critical_quantile = quant_cri_RND)
save(log_ratio_unrelatedness_randomization_test_results, file = './data/analyses/log_ratio_unrelatedness_randomization_test_results.rdata')