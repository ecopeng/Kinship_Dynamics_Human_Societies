rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### load the data
################################################################################
load('./data/analyses/log_ratio_unrelatedness_permutation_test_results.rdata')
load('./data/analyses/log_ratio_unrelatedness_randomization_test_results.rdata')
load('./data/analyses/relatedness_permutation_test_results.rdata')
load('./data/analyses/relatedness_randomization_test_results.rdata')

################################################################################
##### plot test results - patterns of kinship dynamics
################################################################################
res_test_kd <- data.frame(sMSD = as.vector(relatedness_permutation_test_results$distr), approach = 'per')
res_test_kd <- rbind(res_test_kd, data.frame(sMSD = as.vector(relatedness_randomization_test_results$distr), approach = 'rnd'))
res_test_kd$approach <- factor(res_test_kd$approach, levels = c('per', 'rnd'))

vlines_kd <- data.frame(xint = c(relatedness_permutation_test_results$observed_quantile,
                                 relatedness_randomization_test_results$critical_quantile,
                                 relatedness_permutation_test_results$critical_quantile),
                        type = c('obs', 'cri_rnd', 'cri_per'), appr = c('obs', 'rnd', 'per'))

gp <- levels(res_test_kd$approach)
co <- get_ggplot2_colors(length(gp) + 1, 0)

p_test_KD_PER_RND <- ggplot(data = res_test_kd, aes(x = sMSD, fill = factor(approach))) +
  geom_histogram(alpha = .15, bins = 30, linewidth = .5, position = 'identity', color = 'black')

dat_his <- ggplot_build(p_test_KD_PER_RND)$data[[1]]
dat_his$count <- ifelse(dat_his$fill == co[1], -dat_his$count, dat_his$count)
dat_his$fill <- ifelse(dat_his$fill == co[1], gp[1], gp[2])
max_ct <- max(dat_his$count)
min_ct <- min(dat_his$count)

col_bar_vln <- c('obs' = co[1],
                 'rnd' = co[2], 'per' = co[3],
                 'cri_rnd' = co[2], 'cri_per' = co[3])

test_KD <- ggplot(dat_his, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count)) +
  geom_rect(aes(fill = factor(fill)), colour = 'black', linewidth = .1) +
  geom_vline(data = vlines_kd, aes(xintercept = xint, color = type), linetype = 'solid', linewidth = 2.5) +
  scale_y_continuous(limits = c(min_ct, max_ct), labels = abs) +
  labs(x = "sMSD", y = 'count') +
  thm +
  scale_x_continuous(breaks = round(seq(min(dat_his$x), max(dat_his$x), length.out = 8)[c(2, 4, 6)], digits = 2)) +
  theme(legend.position = c(.88, .25),
        legend.direction = 'vertical',
        legend.text = element_text(size = 50, margin = margin(0, 0, 0, 10, 'pt')),
        legend.key.size = unit(2, "lines"),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(colour = "transparent", linewidth = 0),
        legend.box.margin = margin(-5, 0, -5, 0),
        legend.spacing.y = unit(.5, 'cm'),
        axis.text.x = element_text(size = 45),
        axis.title.x = element_text(size = 60),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45),
        plot.margin = unit(c(.25, .25, .05, .25), 'cm')) +
  coord_cartesian(xlim = c(min(dat_his$x), max(dat_his$x))) +
  scale_color_manual(values = col_bar_vln) +
  scale_fill_manual(values = col_bar_vln) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = .75, linewidth = 8)))

################################################################################
##### plot test results - sex differences in kinship dynamics
################################################################################
res_test_sex_diff_kd <- data.frame(sMSD = as.vector(log_ratio_unrelatedness_permutation_test_results$distr), approach = 'per')
res_test_sex_diff_kd <- rbind(res_test_sex_diff_kd, data.frame(sMSD = as.vector(log_ratio_unrelatedness_randomization_test_results$distr), approach = 'rnd'))
res_test_sex_diff_kd$approach <- factor(res_test_sex_diff_kd$approach, levels = c('per', 'rnd'))

vlines_sex_diff_kd <- data.frame(xint = c(log_ratio_unrelatedness_permutation_test_results$observed_quantile,
                                          log_ratio_unrelatedness_randomization_test_results$critical_quantile,
                                          log_ratio_unrelatedness_permutation_test_results$critical_quantile),
                                 type = c('obs', 'cri_rnd', 'cri_per'), appr = c('obs', 'rnd', 'per'))

gp <- levels(res_test_sex_diff_kd$approach)
co <- get_ggplot2_colors(length(gp) + 1, 0)

p_test_DIFF_KD_PER_RND <- ggplot() + 
  geom_histogram(data = res_test_sex_diff_kd, aes(x = sMSD, fill = factor(approach)), alpha = .25, bins = 30, linewidth = .5, position = 'identity', color = 'black')

dat_his <- ggplot_build(p_test_DIFF_KD_PER_RND)$data[[1]]
dat_his$count <- ifelse(dat_his$fill == co[1], -dat_his$count, dat_his$count)
dat_his$fill <- ifelse(dat_his$fill == co[1], gp[1], gp[2])
max_ct <- max(dat_his$count)
min_ct <- min(dat_his$count)

col_bar_vln <- c('obs' = co[1],
                 'rnd' = co[2], 'per' = co[3],
                 'cri_rnd' = co[2], 'cri_per' = co[3])

test_UR <- ggplot(dat_his, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count)) +
  geom_rect(aes(fill = factor(fill)), colour = 'black', linewidth = .1) +
  geom_vline(data = vlines_sex_diff_kd, aes(xintercept = xint, color = type), linetype = 'solid', linewidth = 2.5) +
  scale_y_continuous(limits = c(min_ct, max_ct), labels = abs) +
  labs(x = "sMSD", y = 'count') +
  thm +
  scale_x_continuous(breaks = round(seq(min(dat_his$x), max(dat_his$x), length.out = 8)[c(2, 4, 6)], digits = 2)) +
  theme(legend.position = c(.12, .75),
        legend.direction = 'vertical',
        legend.text = element_text(size = 50, margin = margin(0, 0, 0, 10, 'pt')),
        legend.key.size = unit(2.5, "lines"),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.box.background = element_rect(colour = "transparent", linewidth = 0),
        legend.box.margin = margin(-5, 0, -5, 0),
        legend.spacing.y = unit(.5, 'cm'),
        axis.text.x = element_text(size = 45),
        axis.title.x = element_text(size = 60),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 45),
        plot.margin = unit(c(.25, .25, .05, .25), 'cm')) +
  coord_cartesian(xlim = c(min(dat_his$x), max(dat_his$x))) +
  scale_color_manual(values = col_bar_vln) +
  scale_fill_manual(values = col_bar_vln) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = .75, linewidth = 8)))

################################################################################
##### combine the plots
################################################################################
test_KD <- test_KD + theme(plot.margin = unit(c(.25, .25, .25, 1), "cm"),
                           axis.title.x = element_blank())
test_UR <- test_UR + theme(plot.margin = unit(c(.25, .25, .25, 1), "cm"),
                           legend.position = 'none')

test_res <- ggarrange(test_KD, test_UR, labels = c("(A)", "(B)"),
                      label.args = list(gp = grid::gpar(fontsize = 35, fontface = 2, cex = 2), hjust = .15, vjust = 1.5), ncol = 1)
ggsave(filename = './plots/TEST.pdf', plot = test_res, width = 420, height = 470, units = 'mm')