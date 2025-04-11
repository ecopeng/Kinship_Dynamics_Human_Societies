rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### get the demographic data ready for the model of kinship dynamics
################################################################################
# https://population.un.org/wpp/downloads?folder=Standard%20Projections&group=Mortality
## Number of survivors by age for a hypothetical cohort of 100,000 female newborns
### who would be subject during all their lives to the mortality rates of a given year.
### Data are presented in thousands.
sur_f <- read_excel('./demographics/WPP2024_MORT_F04_3_LIFE_TABLE_SURVIVORS_FEMALE.xlsx', sheet = 1)
sur_f <- as.data.frame(sur_f)
colnames(sur_f) <- sur_f[12, ]
sur_f <- sur_f[13:86, ]
write.csv(sur_f, file = './demographics/sur_f.csv', row.names = F)

### Number of survivors by age for a hypothetical cohort of 100,000 male newborns
### who would be subject during all their lives to the mortality rates of a given year.
### Data are presented in thousands.
sur_m <- read_excel('./demographics/WPP2024_MORT_F04_2_LIFE_TABLE_SURVIVORS_MALE.xlsx', sheet = 1)
sur_m <- as.data.frame(sur_m)
colnames(sur_m) <- sur_m[12, ]
sur_m <- sur_m[13:86, ]
write.csv(sur_m, file = './demographics/sur_m.csv', row.names = F)

### Number of births to women in a particular single age, divided by the number of women in that age.
fer_f <- read_excel('./demographics/WPP2024_FERT_F01_FERTILITY_RATES_BY_SINGLE_AGE_OF_MOTHER.xlsx', sheet = 1)
fer_f <- as.data.frame(fer_f)
colnames(fer_f) <- fer_f[12, ]
fer_f <- fer_f[13:86, ]
write.csv(fer_f, file = './demographics/fer_f.csv', row.names = F)

sur_f <- read.csv('./demographics/sur_f.csv', stringsAsFactors = F)
sur_m <- read.csv('./demographics/sur_m.csv', stringsAsFactors = F)
fer_f <- read.csv('./demographics/fer_f.csv', stringsAsFactors = F)

##### sur_f: [1950, 1986, 2023]
sur_f[, which(colnames(sur_f) == 'X0'):ncol(sur_f)] <- sur_f[, which(colnames(sur_f) == 'X0'):ncol(sur_f)] / 100000
min(sur_f$Year)
round(median(sur_f$Year))
max(sur_f$Year)

##### sur_m: [1950, 1986, 2023]
sur_m[, which(colnames(sur_m) == 'X0'):ncol(sur_m)] <- sur_m[, which(colnames(sur_m) == 'X0'):ncol(sur_m)] / 100000
min(sur_m$Year)
round(median(sur_m$Year))
max(sur_m$Year)

##### fer_f: [1950, 1986, 2023]
fer_f[, which(colnames(fer_f) == 'X15'):ncol(fer_f)] <- fer_f[, which(colnames(fer_f) == 'X15'):ncol(fer_f)] / 1000
min(fer_f$Year)
round(median(fer_f$Year))
max(fer_f$Year)

focal_year <- 1986

sur_f_focal_year <- data.frame(age = as.numeric(gsub('\\D', '', colnames(sur_f[which(sur_f$Year == focal_year), which(colnames(sur_f) == 'X0'):ncol(sur_f)]))),
                               sex = 'female', lx = t(sur_f[which(sur_f$Year == focal_year), which(colnames(sur_f) == 'X0'):ncol(sur_f)][1, ])[, 1])
sur_f_focal_year$qx <- lxtoqx(sur_f_focal_year$lx)
sur_f_focal_year$px <- 1 - sur_f_focal_year$qx

sur_m_focal_year <- data.frame(age = as.numeric(gsub('\\D', '', colnames(sur_m[which(sur_m$Year == focal_year), which(colnames(sur_m) == 'X0'):ncol(sur_m)]))),
                               sex = 'male', lx = t(sur_m[which(sur_m$Year == focal_year), which(colnames(sur_m) == 'X0'):ncol(sur_m)][1, ])[, 1])
sur_m_focal_year$qx <- lxtoqx(sur_m_focal_year$lx)
sur_m_focal_year$px <- 1 - sur_m_focal_year$qx

# reproductive life history stage
fer_f_focal_year <- data.frame(age = as.numeric(gsub('\\D', '', colnames(fer_f[which(fer_f$Year == focal_year), which(colnames(fer_f) == 'X15'):ncol(fer_f)]))),
                               sex = 'female', fx = t(fer_f[which(fer_f$Year == focal_year), which(colnames(fer_f) == 'X15'):ncol(fer_f)][1, ])[, 1])
# pre- & post-reproductive life history stage
fer_f_focal_year <- rbind(fer_f_focal_year, data.frame(age = c(0:(min(fer_f_focal_year$age) - 1), (max(fer_f_focal_year$age) + 1):max(sur_f_focal_year$age)),
                                                       sex = 'female', fx = 0))
fer_f_focal_year <- arrange(fer_f_focal_year, age)

##### write to file for relatedness calculations
fp <- file.path(getwd(), './model_hybrid/par/')
unlink(fp, T)
dir.create(fp)

write.table(sur_f_focal_year$px, file = './model_hybrid/par/SUR_F', sep = '\n', col.names = F, row.names = F)
write.table(sur_f_focal_year$lx, file = './model_hybrid/par/lx_F', sep = '\n', col.names = F, row.names = F)
write.table(sur_m_focal_year$px, file = './model_hybrid/par/SUR_M', sep = '\n', col.names = F, row.names = F)
write.table(sur_m_focal_year$lx, file = './model_hybrid/par/lx_M', sep = '\n', col.names = F, row.names = F)

##### assume males and females have the same fertility curve
write.table(fer_f_focal_year$fx, file = './model_hybrid/par/FEC_F', sep = '\n', col.names = F, row.names = F)
write.table(fer_f_focal_year$fx, file = './model_hybrid/par/FEC_M', sep = '\n', col.names = F, row.names = F)

##### write the numbers and the dispersal rates of females & males, and the rate of local mating to files
# numbers of females and males in each group
NfNm <- expand.grid(N_female = 2:5, N_male = 2:5)
# PAT: patrilocal societies - df, dm, lm
PAT <- c(1, 0, 1)
# MAT: matrilocal societies - df, dm, lm
MAT <- c(0, 1, 1)
# BIL: bilocal societies - df, dm, lm
BIL <- c(.5, .5, 1)
# DUO: duolocal societies - df, dm, lm
DUO <- c(0, 0, 0)

DM <- list(PAT, MAT, BIL, DUO)
count <- 1
for(u in 1:nrow(NfNm)) {
  for(v in DM) {
    write.table(c(NfNm[u, 1], NfNm[u, 2], v), paste0('./model_hybrid/par/NDM_', count - 1),
                sep = '\n', row.names = F, col.names = F)
    cat(NfNm[u, 1], '_', NfNm[u, 2], '_____', v[1], '_', v[2], '_', v[3], '_', '_____', count, '\n', sep = '')
    count <- count + 1
  }
}