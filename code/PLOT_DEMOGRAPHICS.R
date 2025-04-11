rm(list = ls())
arg <- commandArgs(trailingOnly = T)
setwd(arg[1])
source('./code/s.R')

################################################################################
##### load the data
################################################################################
sur_f <- read.table('./model_hybrid/par/SUR_F', sep = '\n', header = F)
colnames(sur_f) <- 'px'
sur_f$sex <- 'F'
sur_f$age <- 1:nrow(sur_f)

sur_m <- read.table('./model_hybrid/par/SUR_M', sep = '\n', header = F)
colnames(sur_m) <- 'px'
sur_m$sex <- 'M'
sur_m$age <- 1:nrow(sur_m)

sur <- rbind(sur_f, sur_m)

fec <- read.table('./model_hybrid/par/FEC_F', sep = '\n', header = F)
colnames(fec) <- 'mx'
fec$age <- 1:nrow(fec)

survival <- ggplot(data = sur, mapping = aes(x = age, y = px, group = sex, color = sex)) +
  geom_line(linewidth = 1, linetype = 'solid', show.legend = T) +
  geom_point(shape = 21, size = 4, fill = 'white') +
  thm +
  scale_y_continuous(breaks = c(0, .4, .8)) +
  scale_x_continuous(breaks = seq(1, 101, by = 50)) +
  labs(x = 'Age', y = 'Survival rate') + 
  theme(legend.position = c(.2, .2),
        legend.direction = 'vertical',
        legend.key = element_rect(fill = "transparent", colour = NA, linewidth = 5),
        legend.text = element_text(size = 30, margin = margin(0, 10, 0, 0, 'pt')),
        legend.box.background = element_rect(colour = "transparent", linewidth = 0),
        legend.box.margin = margin(-5, 0, -5, 0),
        legend.key.size = unit(5, 'line'),
        axis.text.x = element_text(size = 40),
        axis.title.x = element_text(size = 60),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 40),
        plot.margin = unit(c(.25, .25, .05, .25), 'cm'),
        aspect.ratio = .75) +
  guides(color = guide_legend(keywidth = .9, keyheight = .5, default.unit = 'inch', override.aes = list(size = 5))) +
  scale_colour_manual(values = magma(10, begin = .15, end = .85)[c(7, 2)]) +
  scale_fill_manual(values = magma(10, begin = .15, end = .85)[c(7, 2)])

fecundity <- ggplot(data = fec, aes(x = age, y = mx)) +
  geom_line(linewidth = 1) +
  geom_point(shape = 21, size = 4, fill = 'white') +
  thm +
  scale_x_continuous(breaks = seq(1, 101, by = 50)) +
  scale_y_continuous(breaks = c(0, .1, .2)) +
  labs(x = 'Age', y = 'Fertility') +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 40),
        axis.title.x = element_text(size = 60),
        axis.title.y = element_text(size = 60),
        axis.text.y = element_text(size = 40),
        plot.margin = unit(c(.25, .25, .05, .25), 'cm'),
        aspect.ratio = .75)

################################################################################
##### combine the plots
################################################################################
DEMO <- ggarrange(survival, fecundity, labels = c("(A)", "(B)"),
                  label.args = list(gp = grid::gpar(fontsize = 35, fontface = 2, cex = 2), hjust = .5, vjust = 1.5), ncol = 1)
ggsave(filename = './plots/DEMO.pdf', plot = DEMO, width = 360, height = 470, units = 'mm')