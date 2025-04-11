rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

fp <- file.path(getwd(), './plots/')
unlink(fp, T)
dir.create(fp)

################################################################################
##### predicted patterns of kinship dynamics
################################################################################
load('./data/relatedness/PREDICTIONS/relatedness_best_fit.rdata')
r <- relatedness_best_fit[which(relatedness_best_fit$age <= 75), ]
r$soc_N <- paste(r$soc, '(', r$Nf, 'F+', r$Nm, 'M)', sep = '')

p_KD_MODEL <- ggplot(data = r, aes(x = age, y = r_avg, group = interaction(sex, soc), color = sex, fill = sex)) +
  geom_ribbon(aes(ymin = r_avg - sd, ymax = r_avg + sd), alpha = .5, color = 'transparent', show.legend = F) +
  geom_line(linewidth = 3, show.legend = T) +
  thm +
  labs(x = 'Age', y = 'Mean relatedness to group') +
  scale_y_continuous(breaks = function(x) {round(seq(min(x) + (max(x) - min(x)) * .05, max(x), length.out = 30)[c(5, 25)], 2)}) +
  facet_wrap(vars(soc_N), scales = 'free_y', ncol = 2) +
  scale_colour_manual(values = magma(10, begin = .15, end = .85)[c(7, 2)]) +
  scale_fill_manual(values = magma(10, begin = .15, end = .85)[c(7, 2)]) +
  theme(legend.position = c(.75, .1),
        legend.direction = 'vertical',
        legend.key = element_rect(fill = "transparent", colour = NA, linewidth = 5),
        legend.text = element_text(size = 30, margin = margin(0, 10, 0, 0, 'pt')),
        legend.box.background = element_rect(colour = "transparent", linewidth = 0),
        legend.box.margin = margin(-5, 0, -5, 0),
        legend.key.size = unit(3, 'line'),
        strip.background =element_rect(fill = "black", color = "black", linewidth = .75),
        strip.text.x = element_text(size = 30, colour = "white", angle = 0),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 45, margin = margin(0, 5, 0, 0, 'pt'), face = 'bold'),
        axis.text.y = element_text(size = 30),
        plot.margin = unit(c(.25, .25, .05, .25), 'cm'),
        aspect.ratio = .75) +
  guides(color = guide_legend(byrow = T, override.aes = list(linewidth = 4)))

################################################################################
##### predicted patterns of sex differences in kinship dynamics
################################################################################
load('./data/relatedness/PREDICTIONS/log_ratio_unrelatedness_best_fit.rdata')

log_ratio <- log_ratio_unrelatedness_best_fit[which(log_ratio_unrelatedness_best_fit$age <= 75), ]
log_ratio$soc_N <- paste(log_ratio$society, '(', log_ratio$Nf, 'F+', log_ratio$Nm, 'M)', sep = '')
a <- 66

p_log_ratio_ur_f2m_MODEL <- ggplot(log_ratio, aes(x = age, y = log_ur_f2m, color = society)) +
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
  geom_hline(yintercept = 0, linewidth = 1.25, linetype = 'dashed', color = 'black') +
  geom_line(linewidth = 4) +
  thm + 
  scale_y_continuous(limits = function(x){c(min(x), max(x))}, breaks = function(x) {round(seq(min(x) + (max(x) - min(x)) * .05, max(x), length.out = 30)[c(5, 27)], 3)}) + 
  facet_wrap(vars(soc_N), ncol = 2) +
  scale_colour_manual(values = get_ggplot2_colors(length(unique(log_ratio$society)), 0)) +
  labs(x = 'Age', y = 'Sex diff. in mean relatedness to group') +
  theme(legend.position = c(.75, .1),
        legend.direction = 'vertical',
        legend.key = element_rect(fill = "transparent", colour = NA, linewidth = 5),
        legend.text = element_text(size = 30, margin = margin(0, 10, 0, 0, 'pt')),
        legend.box.background = element_rect(colour = "transparent", linewidth = 0),
        legend.box.margin = margin(-5, 0, -5, 0),
        legend.key.size = unit(3, 'line'),
        strip.background =element_rect(fill = "black", color = "black", linewidth = .75),
        strip.text.x = element_text(size = 30, colour = "white", angle = 0),
        axis.text.x = element_text(size = 30),
        axis.title.x = element_text(size = 50),
        axis.title.y = element_text(size = 45, margin = margin(0, 5, 0, 0, 'pt'), face = 'bold'),
        axis.text.y = element_text(size = 30),
        plot.margin = unit(c(.25, .25, .05, .25), 'cm'),
        aspect.ratio = .75) +
  guides(color = guide_legend(byrow = T, override.aes = list(size = 9, linewidth = 4)))

################################################################################
##### combine the plots
################################################################################
p_KD_MODEL <- p_KD_MODEL + theme(plot.margin = unit(c(.75, .25, .25, .75), "cm"))
p_log_ratio_ur_f2m_MODEL <- p_log_ratio_ur_f2m_MODEL + theme(plot.margin = unit(c(.75, .25, .25, .75), "cm"))

p_PRED <- ggarrange(p_KD_MODEL, ggplot() + theme_void(), p_log_ratio_ur_f2m_MODEL, labels = c("(A)", "", "(B)"),
                    label.args = list(gp = grid::gpar(fontsize = 30, fontface = 2, cex = 1.5), hjust = -1, vjust = .5), ncol = 3, widths = c(1, .03, 1))
ggsave(filename = './plots/PRED.pdf', plot = p_PRED, width = 660, height = 340, units = 'mm', device = cairo_pdf)