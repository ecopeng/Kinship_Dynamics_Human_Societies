library(dplyr)
library(stringr)
library(tidyverse)
library(viridis)
library(ggplot2)
library(scales)
library(ggh4x)
library(reshape)
library(proxy)
library(ggarchery)
library(egg)
library(readxl)
library(AMR)
library(patchwork)
library(ggridges)
library(Cairo)
library(tidyr)
library(fmsb)
library(randtests)
library(RColorBrewer)

get_ggplot2_colors <- function(n, prt) {
  if(prt){show_col(hue_pal()(n))}
  return(hue_pal()(n))
}

fun_limits <- function(x){return(c(floor(min(x)), ceiling(max(x))))}
fun_breaks <- function(x){return(seq(floor(min(x)), ceiling(max(x)), by = 1)[c(1, ceiling(max(x)) - floor(min(x)) + 1)])}

thm <- theme(legend.position = "top",
             legend.title = element_blank(),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.direction = "horizontal",
             legend.key = element_rect(fill = "transparent", colour = NA, linewidth = 1),
             legend.text = element_text(size = 17, margin = margin(0, 8, 0, 0, 'pt')),
             legend.key.size = unit(2.5, "lines"),
             legend.margin = margin(0, 0, 0, 0),
             legend.box.margin = margin(-20, 0, -15, 0),
             legend.spacing.x = unit(.5, 'cm'),
             legend.spacing.y = unit(0, 'cm'),
             strip.background = element_rect(fill = "grey0", colour = NA),
             strip.text.x = element_text(color = "white", size = 13, margin = margin(2, 2, 2, 2, "pt"), face = "bold"),
             strip.text.y = element_text(color = "white", size = 13, margin = margin(2, 2, 2, 2, "pt"), face = "bold"),
             panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.background = element_rect(fill = "transparent", colour = NA),
             panel.spacing = unit(.75, "lines"),
             plot.background = element_rect(fill = "transparent", colour = NA),
             axis.line.x = element_line(linewidth = 0, color = NA),
             axis.line.y = element_line(linewidth = 0, color = NA),
             axis.title.x = element_text(size = 15, margin = margin(0, 0, 0, 0, "pt"), face = "bold"),
             axis.title.y = element_text(size = 15, margin = margin(0, 0, 0, 0, "pt"), face = "bold"),
             axis.text.y = element_text(color = "black", size = 12, margin = margin(0, 2, 0, 2, "pt"), face = "bold"),
             axis.text.x = element_text(color = "black", size = 12, margin = margin(3, 0, 3, 0, "pt"), face = "bold"),
             axis.ticks.length = unit(.15, "cm"),
             axis.ticks = element_line(color = "black", linewidth = (.75)),
             aspect.ratio = .618,
             plot.margin = unit(c(-1, .5, -2, 5), "lines"),
             plot.title = element_text(hjust = .5, size = 20))