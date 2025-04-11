rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### load the predicted relatedness data
################################################################################
PAT <- c(1, 0, 1)
MAT <- c(0, 1, 1)
BIL <- c(.5, .5, 1)
DUO <- c(0, 0, 0)
soc <- rbind(PAT, MAT, BIL, DUO)
fn <- list.files('./model_hybrid/out/')
options(dplyr.summarise.inform = F)
r <- data.frame()
for(u in fn) {
  temp <- read.table(paste('./model_hybrid/out/', u, sep = ''), header = T, stringsAsFactors = F)
  ss <- as.numeric(str_extract_all(u, '\\d+([.]\\d+)?')[[1]])
  temp$rep <- ss[1]
  temp$Nf <- ss[2]
  temp$Nm <- ss[3]
  temp$soc <- names(which(apply(soc, 1, FUN = function(x) {identical(x, ss[4:6])})))
  summarised <- temp %>% group_by(sex, age, Nf, Nm, soc) %>%
    dplyr::summarize(r_mean = mean(r2FM), sd = sd(r2FM), n = n(), se = sd(r2FM) / sqrt(n()))
  r <- rbind(r, summarised)
  cat('processing chain [', ss[1], ']\n', sep = '')
}
r$gs <- r$Nf + r$Nm

fp <- file.path(getwd(), './data/relatedness/PREDICTIONS/')
unlink(fp, T)
dir.create(fp)
save(r, file = './data/relatedness/PREDICTIONS/r_PRED.rdata')

################################################################################
##### calcualte log ratio of unrelatedness (female to male)
################################################################################
load('./data/relatedness/PREDICTIONS/r_PRED.rdata')
log_ratio_UR_f2m_PRED <- data.frame()
for(u in unique(r$soc)) {
  for(v in unique(r$Nm)) {
    for(w in unique(r$Nf)) {
      temp <- r[which(r$soc == u & r$Nm == v & r$Nf == w), ]
      temp <- temp %>% group_by(sex, age, Nf, Nm, soc) %>%
        dplyr::summarize(r = mean(r_mean), sd = sd(r_mean), n = n(), se = sd(r_mean) / sqrt(n()))
      a <- sort(unique(temp$age))
      tf <- temp[which(temp$sex == 'F'), ]
      tm <- temp[which(temp$sex == 'M'), ]
      log_ratio_UR_f2m_PRED <- rbind(log_ratio_UR_f2m_PRED,
                                     data.frame(age = unique(temp$age), society = u, Nm = v, Nf = w,
                                                log_ur_f2m = log(
                                                  (1 - tf[match(a, tf$age), ]$r) /
                                                  (1 - tm[match(a, tm$age), ]$r))))
    }
  }
}
save(log_ratio_UR_f2m_PRED, file = './data/relatedness/PREDICTIONS/log_ratio_UR_f2m_PRED.rdata')