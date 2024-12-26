##CREATE FIGURE 2##

library(sf)
library(dplyr)
library(data.table)
library(plotrix)
library(ggplot2)
library(stringr)
library(forcats)
library(ggsci) #Colors for nature journals
library(scales) #To plot palette colors
library(ggh4x) #to group labels by categories https://stackoverflow.com/a/75301076
library(patchwork)
library(Cairo)

rm(list=ls())
dev.off()
gc()


input_dir <- "4A_ANNUAL_STACKS_SUMMARY_DF/ENSEMBLE_SUMMARY"
output_dir <- "FIGURES/FINAL"

grassland <- st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")
grassland_dis <- grassland %>% group_by(ECO_NAME, ECO_ID, G200_REGIO) %>% summarize()
grassland_dis$ID <- seq(1:nrow(grassland_dis))

variables <- dir(input_dir)

axis.df <- data.frame("variables" = variables)
axis.df$title <- c(sort(c(paste("Biovar", 1:19))), variables[20:25])
axis.df$axis <- c(
  "Mean Annual Temperature (\u00b0C)",
  "Mean Temperature of Warmest Quarter (\u00b0C)",
  "Mean Temperature of Coldest Quarter (\u00b0C)",
  "Annual Precipitation (mm)",
  "Precipitation of Wettest Month (mm)",
  "Precipitation of Driest Month (mm)",
  "Precipitation Seasonality (CV)",
  "Precipitation of Wettest Quarter (mm)",
  "Precipitation of Driest Quarter (mm)",
  "Precipitation of Warmest Quarter (mm)",
  "Precipitation of Coldest Quarter (mm)",
  "Annual Mean Diurnal Range (\u00b0C)",
  "Isothermality (%)",
  "Temperature Seasonality (SD*100)",
  "Max Temperature of Warmest Month (\u00b0C)",
  "Min Temperature of Coldest Month (\u00b0C)",
  "Annual Temperature Range (\u00b0C)",
  "Mean Temperature of Wettest Quarter (\u00b0C)",
  "Mean Temperature of Driest Quarter (\u00b0C)",
  "Annual Mean hurs (%)",
  "Annual Precipitation (mm)",
  "Annual Mean rsds (W/m\u00b2)",
  "Annual Mean Temperature (\u00b0C)",
  "Max Temperature of Warmest Month (\u00b0C)",
  "Min Temperature of Coldest Month (\u00b0C)"
)

lp <- list()
lpd <- list()

for (var in 1:length(variables)) {
  
  df.h <- read.csv(list.files(file.path(input_dir, variables[var], "historical"), recursive = TRUE, full.names = TRUE))
  df.h[nrow(df.h)+1,] <- c(nrow(df.h)+1,
                           nrow(df.h)+1,
                           apply(df.h[,3:8], 2, mean),
                           mean(apply(df.h[,3:8], 2, mean)),
                           sd(apply(df.h[,3:8], 2, mean)),
                           std.error(apply(df.h[,3:8], 2, mean)))

  df.585 <- read.csv(list.files(file.path(input_dir, variables[var], "ssp585"), full.names = TRUE))
  df.585[nrow(df.585)+1,] <- c(nrow(df.585)+1,
                           nrow(df.585)+1,
                           apply(df.585[,3:8], 2, mean),
                           mean(apply(df.585[,3:8], 2, mean)),
                           sd(apply(df.585[,3:8], 2, mean)),
                           std.error(apply(df.585[,3:8], 2, mean)))
  
  df.diff <- read.csv(list.files(file.path(input_dir, variables[var], "diff"), full.names = TRUE))
  df.diff[nrow(df.diff)+1,] <- c(nrow(df.diff)+1,
                               nrow(df.diff)+1,
                               apply(df.diff[,3:8], 2, mean),
                               mean(apply(df.diff[,3:8], 2, mean)),
                               sd(apply(df.diff[,3:8], 2, mean)),
                               std.error(apply(df.diff[,3:8], 2, mean)))
  

  df.h$experiment <- "historical"
  
  df.585$experiment <- "ssp585"
  
  df.diff$experiment <- "diff"
  
  df.all <- rbind(df.h, df.585, df.diff)
  
  df.all <- full_join(df.all, st_drop_geometry(grassland_dis), by = "ID")
  
  df.all$ECO_NAME <- ifelse(is.na(df.all$ECO_NAME), "TAE", df.all$ECO_NAME)
  
  df.all$color <- ifelse(df.all$mean < 0, 'negative','positive')
  

  df.all$biorealm <- ifelse(df.all$ECO_NAME %in% c("Santa Marta páramo", "Northern Andean páramo",
                                                   "Cordillera de Merida páramo", "Cordillera Central páramo"), "Neotropic",
                            ifelse(df.all$ECO_NAME %in% c("East African montane moorlands", "Ethiopian montane moorlands"), "Afrotropic",
                                   ifelse(df.all$ECO_NAME == "Central Range sub-alpine grasslands", "Australasia", "")))
  
  df.all$eco_acr <- ifelse(df.all$ECO_NAME == "Santa Marta páramo", "SMP",
                           ifelse(df.all$ECO_NAME == "Northern Andean páramo", "NAP",
                                  ifelse(df.all$ECO_NAME == "Cordillera de Merida páramo", "CMP",
                                         ifelse(df.all$ECO_NAME == "Cordillera Central páramo", "CCP",
                                                ifelse(df.all$ECO_NAME == "East African montane moorlands", "EAMM",
                                                       ifelse(df.all$ECO_NAME == "Ethiopian montane moorlands", "EMM",
                                                              ifelse(df.all$ECO_NAME == "Central Range sub-alpine grasslands", "CRSG",
                                                                     "TAE")))))))
  df.all <- df.all %>%
    mutate(biorealm = fct_relevel(biorealm, "Neotropic", "Afrotropic", "Australasia", ""))
  
  df.all <- df.all %>%
    mutate(ECO_NAME = fct_relevel(ECO_NAME, "Santa Marta páramo", "Northern Andean páramo",
                                  "Cordillera de Merida páramo", "Cordillera Central páramo",
                                  "East African montane moorlands", "Ethiopian montane moorlands",
                                  "Central Range sub-alpine grasslands", "TAE"))
  
  df.all <- df.all %>%
    mutate(eco_acr = fct_relevel(eco_acr, "SMP", "NAP", "CMP", "CCP", "EAMM", "EMM", "CRSG", "TAE"))
  
  #Barplot
  var.plot <- ggplot(df.all[df.all$experiment == "historical" | df.all$experiment == "ssp585",],
                     aes(fill = experiment, y = mean, x = interaction(eco_acr, biorealm, sep = "!"))) + 
    geom_bar(position = "dodge", stat = "identity", width = 0.8, colour = "black") +
    geom_errorbar(aes(ymin = mean-standard_err, ymax = mean+standard_err), width = 0.25, size = 0.35, position = position_dodge(.85)) +
    geom_vline(xintercept = 7.5, linetype = "dashed", size = 0.4) +
    xlab(NULL) + 
    ylab(axis.df$axis[var]) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(guide = guide_axis_nested(delim = "!", extend = 0.75)) +
    scale_fill_manual(values = c("#0072B2BF", "#E69F00BF"),
                      labels = c("1985-2014", "2071-2100")) +
    guides(fill = guide_legend(title = NULL,
                               title.position = "top",
                               title.hjust = 0.5,
                               title.vjust = 0)) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size = 6),
          axis.text = element_text(colour = "black", size = 5),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(size = 8, face = "bold"),
          ggh4x.axis.nesttext.x = element_text(size = 5),
          ggh4x.axis.nestline.x = element_line(colour = "black"))
  
  var.plot
  
  lp[[var]] <- var.plot
  
  if (any(df.all$mean[df.all$experiment == "diff"] <0)) {
    exp <- expansion(mult = .05)
  } else {
    exp <- expansion(mult = 0)
  }
  
  diff.plot <- ggplot(df.all[df.all$experiment == "diff",],
                      aes(fill = color, y = mean, x = interaction(eco_acr, biorealm, sep = "!"))) + 
    geom_bar(position = "dodge", stat = "identity", width = 0.5, colour = "black") +
    geom_errorbar(aes(ymin = mean-standard_err, ymax = mean+standard_err), width = 0.25, size = 0.35, position = position_dodge(.85)) +
    geom_abline(slope = 0, intercept = 0,  col = "black", lty = 1) +
    geom_vline(xintercept = 7.5, linetype = "dashed", size = 0.4) +
    xlab(NULL) + 
    ylab(paste("\u0394", axis.df$axis[var])) +
    scale_y_continuous(expand = exp) +
    scale_x_discrete(guide = guide_axis_nested(delim = "!", extend = 0.75)) +
    scale_fill_manual(values = c(positive = "#0072B2BF", negative = "#E69F00BF")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.title = element_text(size = 6),
          axis.text = element_text(colour = "black", size = 5.5),
          axis.ticks = element_line(colour = "black"),
          plot.title = element_text(size = 8, face = "bold"),
          ggh4x.axis.nesttext.x = element_text(size = 5),
          ggh4x.axis.nestline.x = element_line(colour = "black"))
  
  diff.plot
  
  lpd[[var]] <- diff.plot
  
  # ggsave(plot = diff.plot,
  #        file.path(output_dir, paste0(variables[var], "_diff.pdf")),
  #        device = cairo_pdf,
  #        dpi = 1000,
  #        width = 5.62,
  #        height = 5.21)
}

names(lp) <- axis.df$variables
names(lpd) <- axis.df$variables


# Figure 2 ----------------------------------------------------------------

plot <- lpd$biovar1 + ylab(expression(Delta * "T" * " (°C)")) +
  lpd$biovar12 + ylab(expression(Delta * "P" * " (mm)")) +
  lpd$biovar5 + ylab(expression(Delta * "T"[mw] * " (°C)")) +
  lpd$biovar14 + ylab(expression(Delta * "P"[d] * " (mm)")) +
  lpd$biovar4 + ylab(expression(Delta * "T"[s] * " (SD*100)")) +
  lpd$biovar15 + ylab(expression(Delta * "P"[s] * " (CV)")) +
  plot_layout(ncol = 2, nrow = 3) +
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(size = 8, face = "bold"))

plot

ggsave(plot = plot, "FIGURES/FINAL/FIGURE2.pdf", device = "pdf",  dpi = 1000, width = 130, height = 180, units = "mm") # width = 5.62*2, height = 5.21*2 + 0.19, units = "in"
ggsave(plot = plot, "FIGURES/FINAL/FIGURE2.svg", device = "svg",  dpi = 1000, width = 130, height = 180, units = "mm") # width = 5.62*2, height = 5.21*2 + 0.19, units = "in"
ggsave(plot = plot, "FIGURES/FINAL/FIGURE2.png", device = "png",  dpi = 1000, width = 130, height = 180, units = "mm") # width = 5.62*2, height = 5.21*2 + 0.19, units = "in"


# Figure S2 ---------------------------------------------------------------

plotS3 <- lp$biovar1 + ylab(expression("T" * " (°C)")) +
  lp$biovar12 + ylab(expression("P" * " (mm)")) +
  lp$biovar5 + ylab(expression("T"[mw] * " (°C)")) +
  lp$biovar14 + ylab(expression("P"[d] * " (mm)")) +
  lp$biovar4 + ylab(expression("T"[s] * " (SD*100)")) +
  lp$biovar15 + ylab(expression("P"[s] * " (CV)")) +
  plot_layout(ncol = 2, nrow = 3, guides = "collect") &
  plot_annotation(tag_levels = c("A")) &
  theme(plot.tag = element_text(size = 8, face = "bold"),
        legend.position = "bottom",
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "cm"))


plotS3

ggsave(plot = plotS3, "FIGURES/FINAL/FIGURE_S2.pdf", device = "pdf",  dpi = 1000, width = 130, height = 180, units = "mm")
ggsave(plot = plotS3, "FIGURES/FINAL/FIGURE_S2.svg", device = "svg",  dpi = 1000, width = 130, height = 180, units = "mm")
ggsave(plot = plotS3, "FIGURES/FINAL/FIGURE_S2.png", device = "png",  dpi = 1000, width = 130, height = 180, units = "mm")
