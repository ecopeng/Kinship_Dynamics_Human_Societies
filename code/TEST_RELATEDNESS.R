rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### load the predicted and observed relatedness
################################################################################
load('./data/relatedness/PREDICTIONS/r_PRED.rdata')
r2FM <- r %>% group_by(sex, age, Nf, Nm, soc) %>%
  dplyr::summarize(r_avg = mean(r_mean), sd = sd(r_mean), n = n(), se = sd(r_mean) / sqrt(n()))
r2FM$soc <- factor(r2FM$soc, levels = c('BIL', 'DUO', 'MAT', 'PAT'))

load('./data/relatedness/OBSERVATIONS/r_OBSE.rdata')

################################################################################
##### find the best fit Nm and Nf in the model (sex-age-specific relatedness)
################################################################################
MSD_r <- function(df_DATA, df_PRED) {
  df_DATA <- na.omit(df_DATA)
  idx_in_pred <- match(df_DATA$age, df_PRED$age)
  df_DATA <- df_DATA[!is.na(idx_in_pred), ]
  idx_in_pred <- na.omit(idx_in_pred)
  df_PRED <- df_PRED[idx_in_pred, ]
  DAT <- df_DATA$r_avg
  PRE <- df_PRED$r_avg
  return(sum((PRE - DAT)^2) / length(DAT))
}

v_sum_MSD_r_fm <- function(df_dat, df_mod, society, per_soc_nmnf) {
  vmsd <- c()
  for(p in 1:length(per_soc_nmnf)) {
    df_O <- df_dat[which(df_dat$society == society[p]), ]
    soc_nm_nf <- strsplit(per_soc_nmnf[p], '\\+')[[1]]
    df_P <- df_mod[which(df_mod$soc == soc_nm_nf[1] &
                         df_mod$Nm == as.numeric(soc_nm_nf[2]) &
                         df_mod$Nf == as.numeric(soc_nm_nf[3])), ]
    sum_MSD_fm <- 0
    for(s in c('F', 'M')) {
      sum_MSD_fm <- sum_MSD_fm + MSD_r(df_O[which(df_O$sex == s), ], df_P[which(df_P$sex == s), ])
    }
    vmsd <- c(vmsd, sum_MSD_fm)
  }
  return(vmsd)
}

Nf_Nm <- r2FM[, c('Nm', 'Nf')]
Nf_Nm <- unique(Nf_Nm)
soc_pred <- unique(r2FM$soc)
t_msd_r_fm <- data.frame()

for(v in soc_pred) {
  for(u in 1:nrow(Nf_Nm)) {
    t_pred <- r2FM[which(r2FM$Nm == as.numeric(Nf_Nm[u, 'Nm']) & r2FM$Nf == as.numeric(Nf_Nm[u, 'Nf']) & r2FM$soc == v), ]
    t_obse <- r_OBSE[which(substr(r_OBSE$society, 4, 6) == v), ]
    # each particular soc of given type
    for(t in unique(t_obse$society)) {
      temp <- t_obse[which(t_obse$society == t), c('age', 'sex', 'society', 'r_avg')]
      temp <- na.omit(temp)
      sum_MSD_fm <- 0
      for(s in c('F', 'M')) {
        sum_MSD_fm <- sum_MSD_fm + MSD_r(temp[which(temp$sex == s), ], t_pred[which(t_pred$sex == s), ])
      }
      t_msd_r_fm <- rbind(t_msd_r_fm, data.frame(soc_pred = v,
                                                 Nm = Nf_Nm[u, 'Nm'],
                                                 Nf = Nf_Nm[u, 'Nf'],
                                                 soc_data = t,
                                                 sum_msd_fm = sum_MSD_fm))
    }
  }
}

summary_t_msd_r_fm <- t_msd_r_fm %>% group_by(soc_pred, soc_data) %>% slice(which.min(sum_msd_fm))
save(summary_t_msd_r_fm, file = './data/relatedness/PREDICTIONS/NmNf_relatedness_best_fit.rdata')

relatedness_best_fit <- data.frame()
for(u in 1:nrow(summary_t_msd_r_fm)) {
  temp <- r2FM[which(r2FM$soc == summary_t_msd_r_fm[u, ]$soc_pred &
                     r2FM$Nm == summary_t_msd_r_fm[u, ]$Nm &
                     r2FM$Nf == summary_t_msd_r_fm[u, ]$Nf), ]
  temp$id_soc_data <- summary_t_msd_r_fm[u, ]$soc_data
  relatedness_best_fit <- rbind(relatedness_best_fit, temp)
}
save(relatedness_best_fit, file = './data/relatedness/PREDICTIONS/relatedness_best_fit.rdata')

################################################################################
##### test model predictions: permutation-based approach
################################################################################
fp <- file.path(getwd(), './data/analyses/')
unlink(fp, T)
dir.create(fp)

relatedness_obse <- r_OBSE[, c('age', 'sex', 'society', 'r_avg')]

society <- levels(relatedness_obse$society)
v_soc_NmNf <- paste(summary_t_msd_r_fm$soc_pred, # soc id
                    summary_t_msd_r_fm$Nm,       # Nm
                    summary_t_msd_r_fm$Nf,       # Nf
                    sep = '+')
if(!identical(substr(society, 4, 6), substr(v_soc_NmNf, 1, 3))) {stop('factor levels for *society* inconsistent!')}

tab_soc_NmNf <- table(v_soc_NmNf)
mat_per_soc_nmnf_soc <- t(unique(permut(v_soc_NmNf)))

deno <- 1
for(u in 1:length(tab_soc_NmNf)) {deno <- deno * factorial(as.numeric(tab_soc_NmNf[u]))}
len <- factorial(sum(tab_soc_NmNf)) / deno

MSDs <- matrix(nrow = as.numeric(sum(tab_soc_NmNf)), ncol = ncol(mat_per_soc_nmnf_soc))
for(u in 1:ncol(mat_per_soc_nmnf_soc)) {
  MSDs[, u] <- v_sum_MSD_r_fm(relatedness_obse, relatedness_best_fit, society, mat_per_soc_nmnf_soc[, u])
}

distr_PER <- data.frame(sMSD = colSums(MSDs))
quant_obs_PER <- sum(v_sum_MSD_r_fm(relatedness_obse, relatedness_best_fit, society, v_soc_NmNf))
quant_cri_PER <- quantile(distr_PER, probs = .05, na.rm = T)
cat('PERMUTATION TEST:', '\n', 'actual quantile: ', quant_obs_PER, '\n', 'critical quantile: ', quant_cri_PER, '\n')

relatedness_permutation_test_results <- list(distr = distr_PER, observed_quantile = quant_obs_PER, critical_quantile = quant_cri_PER)
save(relatedness_permutation_test_results, file = './data/analyses/relatedness_permutation_test_results.rdata')

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
  MSDs[, u] <- v_sum_MSD_r_fm(relatedness_obse, relatedness_best_fit, society, soc_nmnf_rnd)
}

distr_RND <- data.frame(sMSD = colSums(MSDs))
quant_obs_RND <- sum(v_sum_MSD_r_fm(relatedness_obse, relatedness_best_fit, society, v_soc_NmNf))
quant_cri_RND <- quantile(distr_RND, probs = .05, na.rm = T)
cat('RANDOMIZATION TEST:', '\n', 'actual quantile: ', quant_obs_RND, '\n', 'critical quantile: ', quant_cri_RND, '\n')

relatedness_randomization_test_results <- list(distr = distr_RND, observed_quantile = quant_obs_RND, critical_quantile = quant_cri_RND)
save(relatedness_randomization_test_results, file = './data/analyses/relatedness_randomization_test_results.rdata')