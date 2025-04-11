rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### MAT & PAT societies (Koster et al)
################################################################################
# summary table Koster et al
tab <- read.csv('./data/relatedness/Koster_et_al_2019/Tab_1.csv')
tab$p_val <- NA
tab$type <- NA

tab$N_native_male <- NA
tab$N_male <- NA
tab$per_cent_native_male <- NA

tab$N_native_female <- NA
tab$N_female <- NA
tab$per_cent_native_female <- NA

# identify patrilocal and matrilocal societies (assuming sexual maturity at age 15)
for(u in 1:nrow(tab)) {
  temp <- read.csv(paste('./data/relatedness/Koster_et_al_2019/', tab$abbreviation[u], '_comm.csv', sep = ''))
  
  N_native_male <- length(which(temp$born == 1 & temp$sex == 1 & temp$age > 15))
  tab$N_native_male[u] <- N_native_male
  N_male <- length(which(temp$sex == 1 & temp$age > 15))
  tab$N_male[u] <- N_male
  tab$per_cent_native_male[u] <- N_native_male / N_male
  
  N_native_female <- length(which(temp$born == 1 & temp$sex == 2 & temp$age > 15))
  tab$N_native_female[u] <- N_native_female
  N_female <- length(which(temp$sex == 2 & temp$age > 15))
  tab$N_female[u] <- N_female
  tab$per_cent_native_female[u] <- N_native_female / N_female
  
  p_val <- prop.test(c(N_native_male, N_native_female), c(N_male, N_female))$p.value
  tab$p_val[u] <- p_val

  tab$type[u] <- ifelse(N_native_male / N_male > N_native_female / N_female, 'PAT', 'MAT')
  if(p_val < .05) {tab$type[u] <- paste(tab$type[u], ' [*]', sep = '')}
}
write.csv(tab, file = './data/relatedness/Koster_et_al_2019/tab.csv')

# load data from those MAT and PAT societies
tab <- tab[str_detect(tab$type, '[*]'), ]
write.csv(tab, file = './data/relatedness/Koster_et_al_2019/tab_mat_pat.csv')

Koster_MAT_PAT <- data.frame()
for(u in 1:nrow(tab)){
  temp <- read.csv(paste('./data/relatedness/Koster_et_al_2019/', tab$abbreviation[u], '_comm.csv', sep = ''))
  temp <- data.frame(society = tab$abbreviation[u], age = temp$age, sex = ifelse(temp$sex - 1, 'F', 'M'), r = temp$avg_r)
  temp <- na.omit(temp)
  temp$age <- round(temp$age)
  summary_temp <- temp %>% group_by(age, sex, society) %>%
    dplyr::summarize(r_avg = mean(r), sd = sd(r), n = n(), se = sd(r) / sqrt(n()))
  Koster_MAT_PAT <- rbind(Koster_MAT_PAT, summary_temp)
}
Koster_MAT_PAT <- Koster_MAT_PAT %>% mutate(
  society = case_when(
    society == 'PU' ~ 'PU_MAT',
    society == 'NI' ~ 'NI_MAT',
    society == 'MA' ~ 'MA_PAT',
    society == 'LA' ~ 'LA_PAT',
    society == 'GA' ~ 'GA_PAT',
    society == 'MP' ~ 'MP_PAT',
    society == 'AZ' ~ 'AZ_PAT',
    society == 'TP' ~ 'TP_PAT',
    society == 'TA' ~ 'TA_PAT'))
Koster_MAT_PAT$society <- factor(Koster_MAT_PAT$society, levels = c('PU_MAT',
                                                                    'NI_MAT',
                                                                    'MA_PAT',
                                                                    'LA_PAT',
                                                                    'GA_PAT',
                                                                    'MP_PAT',
                                                                    'AZ_PAT',
                                                                    'TP_PAT',
                                                                    'TA_PAT'))

################################################################################
##### BIL society (Dyble)
################################################################################
temp <- read.csv(paste('./data/relatedness/Dyble_et_al_2021/ind_data.csv', sep = ''))
Dyble_BIL <- data.frame(society = 'AG_BIL', age = temp$Age, sex = temp$Sex, r = temp$r_to_group)
Dyble_BIL <- na.omit(Dyble_BIL)
Dyble_BIL$age <- round(Dyble_BIL$age)
Dyble_BIL <- Dyble_BIL %>% group_by(age, sex, society) %>%
  dplyr::summarize(r_avg = mean(r), sd = sd(r), n = n(), se = sd(r) / sqrt(n()))

################################################################################
##### DUO society (Wu et al)
################################################################################
### local polynomial regression to get relatedness estimates
Wu_et_al <- data.frame()
### load data females
f2FM <- read.csv('./data/relatedness/Wu_et_al_2013/Wu_et_al_Fig4a_blue.csv', header = F)
colnames(f2FM) <- c('age', 'r_range')
f2FM$sex <- 'F'
f2FM$age <- round(f2FM$age)
summary_f2FM <- f2FM %>% group_by(age, sex) %>%
  dplyr::summarize(r = mean(r_range))
Wu_et_al <- rbind(Wu_et_al, summary_f2FM)
### load data males
m2FM <- read.csv('./data/relatedness/Wu_et_al_2013/Wu_et_al_Fig4c_blue.csv', header = F)
colnames(m2FM) <- c('age', 'r_range')
m2FM$sex <- 'M'
m2FM$age <- round(m2FM$age)
summary_m2FM <- m2FM %>% group_by(age, sex) %>%
  dplyr::summarize(r = mean(r_range))
Wu_et_al <- rbind(Wu_et_al, summary_m2FM)

### estimates female relatedness
age <- Wu_et_al[which(Wu_et_al$sex == 'F'), ]$age
r <- Wu_et_al[which(Wu_et_al$sex == 'F'), ]$r
l <- loess(r ~ age, span = .75)
age_seq <- seq(min(Wu_et_al$age), max(Wu_et_al$age), by = 1)
f2_pred <- predict(l, age_seq)
fitted_data_f <- data.frame(age = age_seq, r = f2_pred, sex = 'F')

### estimates male relatedness
age <- Wu_et_al[which(Wu_et_al$sex == 'M'), ]$age
r <- Wu_et_al[which(Wu_et_al$sex == 'M'), ]$r
l <- loess(r ~ age, span = .75)
age_seq <- seq(min(Wu_et_al$age), max(Wu_et_al$age), by = 1)
m2_pred <- predict(l, age_seq)
fitted_data_m <- data.frame(age = age_seq, r = m2_pred, sex = 'M')

Wu_DUO <- rbind(fitted_data_f, fitted_data_m)
Wu_DUO$society <- 'MS_DUO'

################################################################################
##### combine all the datasets: BIL(1), DUO(1), MAT(2), PAT(7)
################################################################################
soc <- c('AG_BIL',
         'MS_DUO',
         'PU_MAT', 'NI_MAT',
         'MA_PAT', 'LA_PAT', 'GA_PAT', 'MP_PAT', 'AZ_PAT', 'TP_PAT', 'TA_PAT')
colnames(Wu_DUO)[which(colnames(Wu_DUO) == 'r')] <- 'r_avg'
Wu_DUO$sd <- NA
Wu_DUO$n <- 1
Wu_DUO$se <- NA

r_OBSE <- rbind(Koster_MAT_PAT, Dyble_BIL, Wu_DUO)
r_OBSE$society <- factor(r_OBSE$society, levels = soc)

fp <- file.path(getwd(), './data/relatedness/OBSERVATIONS/')
unlink(fp, T)
dir.create(fp)

save(r_OBSE, file = './data/relatedness/OBSERVATIONS/r_OBSE.rdata')

log_ratio_UR_f2m_OBSE <- data.frame()
for(u in as.character(unique(r_OBSE$society))) {
  t <- r_OBSE[which(r_OBSE$society == u), ]
  a <- sort(unique(t$age))
  tf <- t[which(t$sex == 'F'), ]
  tm <- t[which(t$sex == 'M'), ]
  log_ratio_UR_f2m_OBSE <- rbind(log_ratio_UR_f2m_OBSE,
                                 data.frame(age = a, society = u, log_ur_f2m = log(
                                   (1 - tf[match(a, tf$age), ]$r_avg) / (1 - tm[match(a, tm$age), ]$r_avg))))
}
log_ratio_UR_f2m_OBSE$society <- factor(log_ratio_UR_f2m_OBSE$society, levels = soc)
log_ratio_UR_f2m_OBSE$soc <- substr(log_ratio_UR_f2m_OBSE$society, 4, 6)
log_ratio_UR_f2m_OBSE$soc <- factor(log_ratio_UR_f2m_OBSE$soc, levels = c('BIL', 'DUO', 'MAT', 'PAT'))

save(log_ratio_UR_f2m_OBSE, file = './data/relatedness/OBSERVATIONS/log_ratio_UR_f2m_OBSE.rdata')