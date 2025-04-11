rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### observed patterns of kinship dynamics
################################################################################
load('./data/relatedness/OBSERVATIONS/r_OBSE.rdata')
soc <- c('AG_BIL',
         'MS_DUO',
         'PU_MAT', 'NI_MAT',
         'MA_PAT', 'LA_PAT', 'GA_PAT', 'MP_PAT', 'AZ_PAT', 'TP_PAT', 'TA_PAT')
design <- "A######
           B######
           CD#####
           EFGHIJK"
r_OBSE$society <- factor(r_OBSE$society, levels = soc)

p_KD_DATA <- ggplot() +
  geom_point(r_OBSE, mapping = aes(x = age, y = r_avg, color = sex, fill = sex), size = 2.5, shape = 21, fill = 'transparent', alpha = .7, show.legend = T) +
  geom_smooth(method = loess, r_OBSE[which(r_OBSE$society != 'MS_DUO'), ], mapping = aes(x = age, y = r_avg, color = sex, fill = sex), linewidth = 3, show.legend = F, se = T) +
  geom_smooth(method = loess, r_OBSE[which(r_OBSE$society != 'MS_DUO'), ], mapping = aes(x = age, y = r_avg, color = sex, fill = sex), alpha = 1, show.legend = T, se = F) +
  thm + 
  labs(x = 'Age', y = 'Mean relatedness to group') +
  facet_manual(vars(society), design = design, scales = 'free_y') +
  theme(legend.position = c(.95, .75),
        legend.direction = "vertical",
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 30, margin = margin(0, 0, 0, 0, 'pt')),
        legend.key.size = unit(5, "lines"),
        legend.spacing.y = unit(-1.25, 'cm'),
        strip.text.x = element_text(size = 25, colour = "white", angle = 0),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 50),
        axis.text.y = element_text(size = 30),
        plot.margin = unit(c(.25, .25, .05, .25), 'cm'),
        aspect.ratio = .75) +
  scale_fill_manual(values = magma(10, begin = .15, end = .85)[c(7, 2)]) +
  scale_colour_manual(values = magma(10, begin = .15, end = .85)[c(7, 2)]) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 9, linewidth = 4))) +
  scale_y_continuous(limits = function(x){c(min(x), max(x))}, breaks = function(x) {round(seq(min(x) + (max(x) - min(x)) * .05, max(x), length.out = 30)[c(5, 27)], 3)})
ggsave(filename = './plots/OBSE_KD_DATA.pdf', plot = p_KD_DATA, width = 800, height = 300, units = 'mm', device = cairo_pdf)

################################################################################
##### observed patterns of sex differences in kinship dynamics
################################################################################
load('./data/relatedness/OBSERVATIONS/log_ratio_UR_f2m_OBSE.rdata')
dat <- na.omit(log_ratio_UR_f2m_OBSE)

col <- c(brewer.pal(n = 9, 'Reds')[5], brewer.pal(n = 9, 'YlGn')[5], brewer.pal(n = 9, 'Blues')[c(6, 8)], brewer.pal(n = 9, 'RdPu')[3:9])
show_col(col)
a <- 80

p_log_ratio_ur_f2m_DATA <- ggplot(dat, aes(x = age, y = log_ur_f2m, group = interaction(society, soc), color = society, fill = society)) + 
  geom_segment(x = a, y = .01, xend = a, yend = .07,
               lineend = 'butt', linejoin = 'mitre', size = 2,
               arrow = arrow(length = unit(.15, 'inches'), type = 'closed'), color = 'black', linewidth = 1) +
  annotate('text', label = '\u2642>\u2640', x = a, y = .1, size = 12,
           colour = 'black', fontface = 2, family = 'serif') +
  geom_segment(x = a, y = -.01, xend = a, yend = -.07,
               lineend = 'butt', linejoin = 'mitre', size = 2,
               arrow = arrow(length = unit(.15, 'inches'), type = 'closed'), color = 'black', linewidth = 1) +
  annotate('text', label = '\u2640>\u2642', x = a, y = -.1, size = 12,
           colour = 'black', fontface = 2, family = 'serif') +
  geom_point(size = 5, shape = 21) +
  geom_smooth(method = loess, linewidth = 4, show.legend = T, se = T, alpha = .1) +
  geom_hline(yintercept = 0, linewidth = 1.5, linetype = 'dashed', color = 'black') +
  facet_wrap(vars(soc)) +
  thm + 
  scale_y_continuous(limits = function(x){c(min(x), max(x))}, breaks = function(x) {round(seq(min(x) + (max(x) - min(x)) * .05, max(x), length.out = 30)[c(5, 27)], 3)}) + 
  labs(x = 'Age', y = 'Sex diff. in mean relatedness to group') +
  theme(legend.direction = "horizontal",
        # legend.position = c(.95, .75),
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.text = element_text(size = 30, margin = margin(0, 8, 0, 0, 'pt')),
        legend.key.size = unit(3, "lines"),
        legend.margin = margin(4, 4, 15, 4),
        legend.spacing.x = unit(1, 'cm'),
        strip.background =element_rect(fill = "black", color = "black", linewidth = .75),
        strip.text.x = element_text(size = 30, colour = "white", angle = 0),
        axis.text.x = element_text(size = 35),
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 47),
        axis.text.y = element_text(size = 35),
        plot.margin = unit(c(.25, .25, .25, .25), "cm"),
        aspect.ratio = .75) +
  scale_colour_manual(values = col) +
  scale_fill_manual(values = col) +
  guides(color = guide_legend(byrow = T, nrow = 2, override.aes = list(size = 7, linewidth = 4)),
         fill = guide_legend(nrow = 2, byrow = T))
ggsave(filename = './plots/OBSE_KD_DIFF.pdf', plot = p_log_ratio_ur_f2m_DATA, width = 500, height = 430, units = 'mm', device = cairo_pdf)