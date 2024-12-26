##PRODUCE HYDROCLIMATIC SPACES##

library(sf)
library(dplyr)
library(ggplot2)
library(ggdensity)
library(purrr)
library(forcats)
library(patchwork)
library(stringr)
library(raster)
library(rgdal)
library(gmodels)
library(rnaturalearth)
library(tidyr)
library(ggspatial)
library(Cairo)
library(scales) #To plot palette colors
library(ggthemes)
library(cowplot)
library(ggpubr)
library(measurements)
library(units)
library(stringr)
library(ggsci)
library(ggpattern)
library(ggh4x) #to group labels by categories https://stackoverflow.com/a/75301076

##https://gis.stackexchange.com/questions/444689/create-density-polygons-from-spatial-points-in-r

rm(list = ls())
dev.off()
gc()

input <- "4A_ANNUAL_STACKS_SUMMARY_DF/INTERPOLATED"

# Load and prepare data ---------------------------------------------------

#Modify according to variables being used to create the hydroclimatic space
#Any of long_term, extreme or seasonality!
hspaces <- c("long_term", "extreme", "seasonality")
variables <- dir(input)
models <- dir(file.path(input, "biovar1"))
hs_df <- data.frame(matrix(ncol = 6, nrow = 0))

for (hspace in hspaces) {
  
  mod_df <- data.frame(matrix(ncol = 6, nrow = 0))
  
  if (hspace == "long_term") {
    
    varx <- "biovar1"
    vary <- "biovar12"
    
    labx <- "Temperature (\u00b0C)"
    laby <- "Precipitation (mm)"
    
  } else if (hspace == "extreme") {
    
    varx <- "biovar5"
    vary <- "biovar14"
    
    labx <- "Temperature (\u00b0C)"
    laby <- "Precipitation (mm)"
    
  } else {
    
    varx <- "biovar4"
    vary <- "biovar15"
    
    labx <- "Temperature SD * 100 (\u00b0C)"
    laby <- "Precipitation CV (%)"
  }
  
  #SCENARIO I
  # Santa Marta Páramo is joined to northern andean páramo for its lower number of pixels coverage
  # Load ecoregions .shp and dissolve by ecoregion
  ecoregions <- st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")
  ecoregions[which(ecoregions$ECO_NAME == "Santa Marta páramo"),][4] <- "Northern Andean páramo"
  ecoregions_dis <- ecoregions %>% group_by(ECO_NAME, G200_REGIO) %>% summarize()
  ecoregions_dis$ID <- seq(1:nrow(ecoregions_dis))
  
  for (mod in models) {
    
    # Load and prepare data frame with pixel values for varx and vary and for h and 585
    
    df.h <- full_join(read.csv(file.path(input, varx, mod, "historical", paste0(varx, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "cen_1985-2014.csv")), encoding = "UTF-8"),
                      read.csv(file.path(input, vary, mod, "historical", paste0(vary, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "cen_1985-2014.csv")), encoding = "UTF-8"))
    
    df.h <- df.h[-c(1,2,6)]
    names(df.h) <- c("ecoregion", "G200_REGION", "polygon_id", "ID", "X", "Y")
    df.h$exp <- "historical"
    df.h$ecoregion[which(df.h$ecoregion == "Santa Marta páramo")] <- "Northern Andean páramo"
    
    df.585 <- full_join(read.csv(file.path(input, varx, mod, "SSP585", paste0(varx, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "cen_2071-2100.csv")), encoding = "UTF-8"),
                        read.csv(file.path(input, vary, mod, "SSP585", paste0(vary, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "cen_2071-2100.csv")), encoding = "UTF-8"))
    
    df.585 <- df.585[-c(1,2,6)]
    names(df.585) <- c("ecoregion", "G200_REGION", "polygon_id", "ID", "X", "Y")
    df.585$exp <- "ssp585"
    df.585$ecoregion[which(df.585$ecoregion == "Santa Marta páramo")] <- "Northern Andean páramo"
    
    df <- rbind(df.h, df.585)
    
    df$eco <- ifelse(df$ecoregion == "Central Range sub-alpine grasslands", "CRSG",
                     ifelse(df$ecoregion == "East African montane moorlands", "EAMM",
                            ifelse(df$ecoregion == "Ethiopian montane moorlands", "EMM",
                                   ifelse(df$ecoregion == "Northern Andean páramo", "NAP",
                                          ifelse(df$ecoregion == "Cordillera de Merida páramo", "CMP",
                                                 ifelse(df$ecoregion == "Cordillera Central páramo", "CCP", NA))))))
    
    rm(df.h, df.585)
    
    # Prepare xlim and ylim for kdensity plots
    
    if (min(df$X) < 0) {
      
      xmin <- floor(min(df$X)) -
        floor(min(df$X) * 0.1)
      
    } else {
      
      xmin <- 0
    }
    
    xmax <- ceiling(max(df$X)) +
      ceiling(max(df$X) * 0.1)
    
    xlim <- c(xmin, xmax)
    
    
    if (min(df$Y) < 0) {
      
      ymin <- floor(min(df$Y)) -
        floor(min(df$Y) * 0.1)
      
    } else {
      
      ymin <- 0
    }
    
    ymax <- ceiling(max(df$Y)) +
      ceiling(max(df$Y) * 0.1)
    
    ylim <- c(ymin, ymax)
    
    # By ecoregion - kdensity fig ------------------------------------------------------------
    
    # Create empty variables to store data while looping
    pixels.intersect <- data.frame(row.names = names(df)) # Variable storing all the intersected pixels
    pixels.diff <- data.frame(row.names = names(df)) # Variable storing all the not intersected pixels
    lplots <- list() # List to store the plots
    lpolygons.h <- list() # List to store the historical polygons
    lpoints.585 <- list() # List to store the 585 polygons
    lpoints.intersect <- list() # List to store the intersected points
    lpoints.diff <- list() # List to store the not intersected points
    int <- list() #List to store the movement intensity
    dir <- list() #List to store the movement direction
    sev <- list() #List to store the movement severity
    
    # Loop by ecoregion to produce kdensity plot
    
    for (eco in unique(df$ecoregion)) {
      
      # Get centroid of the points and prepare a data frame for later plotting
      centroid.h <- matrix(colMeans(df[df$exp == "historical" & df$ecoregion == eco,][c("X","Y")]))
      centroid.h <- data.frame(X <- centroid.h[1], Y <- centroid.h[2])
      names(centroid.h) <- c("X", "Y")
      
      centroid.585 <- matrix(colMeans(df[df$exp == "ssp585" & df$ecoregion == eco,][c("X","Y")]))
      centroid.585 <- data.frame(X <- centroid.585[1], Y <- centroid.585[2])
      names(centroid.585) <- c("X", "Y")
      
      centroid <- rbind(centroid.h, centroid.585)
      centroid$exp <- c("historical", "ssp585")
      centroid$group <- 1
      
      # Compute magnitude (intensity) and direction of the vector of change
      
      min_x <- min(df$X)
      max_x <- max(df$X)
      
      min_y <- min(df$Y)
      max_y <- max(df$Y)
      
      V <- c(centroid.585$X, centroid.585$Y) - c(centroid.h$X, centroid.h$Y)
      names(V) <- c("x", "y")
      
      V_norm <- c((centroid.585$X - min_x)/(max_x - min_x),
                  (centroid.585$Y - min_y)/(max_y - min_y)) - 
        c((centroid.h$X - min_x)/(max_x - min_x),
          (centroid.h$Y - min_y)/(max_y - min_y))
      names(V_norm) <- c("x", "y")
      
      intensity <- round(sqrt(sum(V_norm^2)), 2)
      
      direction <- round(atan2(V_norm["y"], V_norm["x"]) * (180/pi), 2)
      
      # Create the basic kdensity plot
      plot <- ggplot() +
        geom_hdr(data = df[df$ecoregion == eco & df$exp == "historical",], aes(X, Y, fill = "#0072B2"), probs = 0.95, method = "kde", alpha = 0.5, color = "#0A0A0ABB", size = 0.2) +
        geom_hdr(data = df[df$ecoregion == eco & df$exp == "ssp585",], aes(X, Y, fill = "#E69F00"), probs = 0.95, method = "kde", alpha = 0.5, color = '#0A0A0ABB', size = 0.2) +
        geom_line(data = centroid, aes(X, Y, group = group), size = 0.5, color = "black", linetype = "solid", arrow = arrow(type = "open", length = unit(0.075, "inches"), angle = 30)) +
        annotate(geom = "text", x = xlim[2]*0.15, y = ylim[2]*0.92, label = paste0("\u03b8: ", direction, "\u00b0"), size = 2.25, fontface ="bold") +
        annotate(geom = "text", x = xlim[2]*0.15, y = ylim[2]*0.85, label = paste0("I: ", intensity), size = 2.25, fontface = "bold") +
        scale_x_continuous(limits = xlim, breaks = breaks_extended(6), expand = c(0,0)) +
        scale_y_continuous(limits = ylim, breaks = breaks_extended(6), expand = c(0,0)) +
        scale_fill_manual(values = c("#0072B2", "#E69F00"),
                          labels = c("1985-2014", "2071-2100")) +
        xlab(labx) + ylab(laby) +
        theme_classic() +
        theme(legend.position = "none",
              axis.text = element_text(colour = "Black", size = 7),
              axis.title = element_text(size = 8),
              axis.line = element_line(colour = 'black', size = 0.1),
              element_line(colour = "black", size = 0.1),
              plot.title = element_text(size = 8))
      
      plot
      
      # Extract the kdensity geometries from the plot and create sf polygon objects
      data2d.h <- ggplot2::layer_data(plot, i = 1)
      polygon.h <- st_as_sf(data2d.h, coords = c("x", "y")) %>%
        group_by(subgroup) %>%
        summarise(geometry = st_combine(geometry)) %>%
        st_cast("POLYGON")
      
      data2d.585 <- ggplot2::layer_data(plot, i = 2)
      polygon.585 <- st_as_sf(data2d.585, coords = c("x", "y")) %>%
        group_by(subgroup) %>%
        summarise(geometry = st_combine(geometry)) %>%
        st_cast("POLYGON")
      
      # Get % of reference hydroclimatic space intersecting future hydroclimatic space
      
      if (any(st_intersects(polygon.h, polygon.585, sparse = FALSE))) {
        
        inters <- st_intersection(polygon.h, polygon.585)
        
        severity <- round(100 - ((sum(st_area(inters)) / sum(st_area(polygon.h))) * 100), 2)
        
      } else {
        
        severity <- 100 
        
      }
      
      int[[eco]] <- intensity
      dir[[eco]] <- direction
      sev[[eco]] <- severity
      
      # Create sf point objects using the x/y space coordinates for 585
      points.585 <- st_as_sf(df[df$ecoregion == eco & df$exp == "ssp585",], coords = c("X", "Y"))
      points.585$X <- df$X[df$ecoregion == eco & df$exp == "ssp585"]
      points.585$Y <- df$Y[df$ecoregion == eco & df$exp == "ssp585"]
      
      # Get the points from 585 which intersect the polygon from h in the x/y space
      points.intersect <- st_filter(points.585, polygon.h)
      
      if (nrow(points.intersect) > 0) {
        
        points.intersect$intersect <- "1"
        
      }
      
      
      # Get the points from 585 which don't intersect the polygon from h in the x/y space
      points.diff <- points.585[!(points.585$ID %in% points.intersect$ID),]
      
      if (nrow(points.diff) > 0) {
        
        points.diff$intersect <- "0"
        
      }
      
      # Store iteratively the points intersected or not by ecoregion
      
      pixels.intersect <- rbind(pixels.intersect, points.intersect)
      
      
      pixels.diff <- rbind(pixels.diff, points.diff)
      
      
      pixels.final <- rbind(points.intersect, points.diff)
      
      pixels.final$intersect <- factor(pixels.final$intersect, levels = c("0", "1"), ordered = TRUE)
      
      
      col <- c("0" = "#E69F00", "1" = "#0072B2")
      
      # Add intersected and not intersected points to the kdensity plot
      plot <- plot + 
        geom_point(data = pixels.final, aes(X, Y, color = intersect), size = 0.1, shape = 16) +
        scale_color_manual(guide = guide_legend(reverse = TRUE), values = col, labels = c("Breached areas", "Safe areas")) +
        annotate(geom = "text", x = xlim[2]*0.15, y = ylim[2]*0.78, label = paste0("S: ", severity, "%"), size = 2.25, fontface = "bold") +
        labs(colour = NULL,
             fill = NULL) +
        theme(legend.position = "none")
      
      plot
      
      if (eco == "Central Range sub-alpine grasslands" | eco == "Northern Andean páramo") {
        
        if (hspace == "long_term") {
          
          plot <- plot + ggtitle("LONG-TERM") + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
          
        } else if (hspace == "extreme") {
          
          plot <- plot + ggtitle("EXTREMES") + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
          
        } else  {
          
          plot <- plot + ggtitle("SEASONALITY") + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
          
        }
        
      } else {
        
        plot <- plot + ggtitle("") + theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5))
        
      }
      
      xdens <- axis_canvas(plot, axis = "x") +
        geom_density(data = df[df$ecoregion == eco,], aes(x = X, fill = exp),
                     alpha = 0.5, size = 0.2) +
        fill_palette(c("#0072B2", "#E69F00"))+
        geom_rug(data = pixels.final, aes(x = X, color = intersect), length = unit(0.1, "cm")) +
        color_palette(col)
      
      ydens <- axis_canvas(plot, axis = "y", coord_flip = TRUE) +
        geom_density(data = df[df$ecoregion == eco,], aes(x = Y, fill = exp),
                     alpha = 0.5, size = 0.2) +
        fill_palette(c("#0072B2", "#E69F00")) +
        geom_rug(data = pixels.final, aes(x = Y, color = intersect), length = unit(0.1, "cm")) +
        color_palette(col) +
        coord_flip()
      ydens
      
      p1 <- insert_xaxis_grob(plot, xdens, grid::unit(.2, "null"), position = "top")
      p2 <- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
      ggdraw(p2)
      
      # Store variables iteratively by ecoregion when model == "ENSEMBLE"
      
      if (mod == "ENSEMBLE") {
        
        lplots[[eco]] <- ggdraw(p2)
        lpolygons.h[[eco]] <- polygon.h
        lpoints.585[[eco]] <- points.585
        lpoints.intersect[[eco]] <- points.intersect
        lpoints.diff[[eco]] <- points.diff
        
      }
      
    }
    
    legend <- get_legend(
      plot + 
        guides(colour = guide_legend(nrow = 1, override.aes = list(size = 1.5))) +
        theme(legend.position = "bottom",
              legend.text = element_text(size = 8),
              legend.key.size = unit(0.5, "cm")))
    
    if (mod == "ENSEMBLE") {
      
      save(lplots, lpolygons.h, lpoints.585, lpoints.intersect, lpoints.diff, legend,
           pixels.diff, pixels.intersect, pixels.final,
           file = paste0("FIGURES/FINAL/", hspace, "_data.RData"))
      
      rm(lplots, lpolygons.h, lpoints.585, lpoints.intersect, lpoints.diff, legend,
         pixels.diff, pixels.intersect, pixels.final)
      
    } else {
      
      mod_df <- rbind(mod_df,
                      cbind(rep(hspace, 6), rep(mod, 6), unique(ecoregions$ECO_NAME),
                            unlist(int, use.names = FALSE), unlist(sev, use.names = FALSE),
                            unlist(dir, use.names = FALSE)))
      
    }
  }
  
  colnames(mod_df) <- c("hspace","model", "ecoregion", "intensity", "severity", "direction")
  
  hs_df <- rbind(hs_df, mod_df)
  
  if (which(hspaces %in% hspace) == length(hspaces)) {
    
    save(hs_df, file = "FIGURES/FINAL/hspaces_mod_data.RData")
    
  }
  
}


# Compute 95% CI ----------------------------------------------------------

# load("FIGURES/FINAL/hspaces_mod_data.RData")

ci_df <- data.frame(matrix(ncol = 7, nrow = 0))

for (hspace in unique(hs_df$hspace)) {
  
  for (eco in unique(hs_df$ecoregion)) {
    
    for (var in names(hs_df[4:6])) {
      
      ci <- gmodels::ci(as.numeric(hs_df[hs_df$hspace == hspace & hs_df$ecoregion == eco, var]), confidence = 0.95)
      
      ci_df <- rbind(ci_df, c(hspace, eco, var, round(as.numeric(ci), 4)))
      
    }
    
  }
  
}

names(ci_df) <- c("hspace", "ecoregion", "variable", "Estimate", "CI lower", "CI upper", "Std. Error")
write.csv(ci_df, file = "FIGURES/FINAL/95_ci.csv", fileEncoding = "UTF-8")
write_xlsx(ci_df, "FIGURES/FINAL/95_ci.xlsx")

# Join the plots to produce figure 3 --------------------------------------
#cowplot

load("FIGURES/FINAL/long_term_data.RData")
lplots.lt <- lplots

load("FIGURES/FINAL/extreme_data.RData")
lplots.ex <- lplots

load("FIGURES/FINAL/seasonality_data.RData")
lplots.se <- lplots


p <- plot_grid(
  lplots.lt$`Central Range sub-alpine grasslands` + theme(legend.position = "none"),
  lplots.ex$`Central Range sub-alpine grasslands` + theme(legend.position = "none"),
  lplots.se$`Central Range sub-alpine grasslands` + theme(legend.position = "none"),
  lplots.lt$`Cordillera de Merida páramo` + theme(legend.position="none"),
  lplots.ex$`Cordillera de Merida páramo` + theme(legend.position="none"),
  lplots.se$`Cordillera de Merida páramo` + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "", "", "B", "", ""),
  label_size = 10,
  nrow = 2,
  ncol = 3
)

p

plegend <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, .1))

plegend

ggsave(plot = plegend, paste0("FIGURES/FINAL/FIGURE3.pdf"), device = "pdf", dpi = 1000, width = 180, height = 120, units = "mm")
ggsave(plot = plegend, paste0("FIGURES/FINAL/FIGURE3.svg"), device = "svg", dpi = 1000, width = 180, height = 120, units = "mm")
ggsave(plot = plegend, paste0("FIGURES/FINAL/FIGURE3.png"), device = "png", dpi = 1000, width = 180, height = 120, units = "mm")



# Join the plots to produce figure S3
#cowplot

p <- plot_grid(
  lplots.lt$`Northern Andean páramo`+ theme(legend.position = "none"),
  lplots.ex$`Northern Andean páramo` + theme(legend.position = "none"),
  lplots.se$`Northern Andean páramo` + theme(legend.position = "none"),
  lplots.lt$`Cordillera Central páramo` + theme(legend.position="none"),
  lplots.ex$`Cordillera Central páramo`  + theme(legend.position="none"),
  lplots.se$`Cordillera Central páramo`  + theme(legend.position="none"),
  lplots.lt$`Ethiopian montane moorlands` + theme(legend.position="none"),
  lplots.ex$`Ethiopian montane moorlands`  + theme(legend.position="none"),
  lplots.se$`Ethiopian montane moorlands`  + theme(legend.position="none"),
  lplots.lt$`East African montane moorlands`+ theme(legend.position="none"),
  lplots.ex$`East African montane moorlands`  + theme(legend.position="none"),
  lplots.se$`East African montane moorlands` + theme(legend.position="none"),
  align = 'vh',
  labels = c("A", "", "", "B", "", "", "C", "", "", "D", "", ""),
  label_size = 10,
  nrow = 4,
  ncol = 3
)

p

plegend <- plot_grid(p, legend, ncol = 1, rel_heights = c(1, .1))

plegend

ggsave(plot = plegend, paste0("FIGURES/FINAL/FIGURE_S3.pdf"), device = "pdf", dpi = 1000, width = 180, height = 240, units = "mm")
ggsave(plot = plegend, paste0("FIGURES/FINAL/FIGURE_S3.svg"), device = "svg", dpi = 1000, width = 180, height = 240, units = "mm")
ggsave(plot = plegend, paste0("FIGURES/FINAL/FIGURE_S3.png"), device = "png", dpi = 1000, width = 180, height = 240, units = "mm")

# By ecoregion - Maps --------------------------------------------------------------------

#Prepare data
ecoregions <- st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")
ecoregions[which(ecoregions$ECO_NAME == "Santa Marta páramo"),][4] <- "Northern Andean páramo"
ecoregions_dis <- ecoregions %>% group_by(ECO_NAME, G200_REGIO) %>% summarize()
ecoregions_dis$ID <- seq(1:nrow(ecoregions_dis))

eco.sam <- subset(ecoregions, ecoregions$ECO_NAME == "Cordillera Central páramo" |
                    ecoregions$ECO_NAME == "Northern Andean páramo" |
                    ecoregions$ECO_NAME == "Cordillera de Merida páramo")

box.sam <- st_bbox(eco.sam)

eco.afr <- subset(ecoregions, ecoregions$ECO_NAME == "East African montane moorlands" |
                    ecoregions$ECO_NAME == "Ethiopian montane moorlands")

box.afr <- st_bbox(eco.afr)

eco.png <- subset(ecoregions, ecoregions$ECO_NAME == "Central Range sub-alpine grasslands")

world <- ne_countries(scale = "medium", returnclass = "sf")


# Prepare polygons to highlight ecoregions --------------------------------

grid <- st_read("GIS/GRID/GRID.shp")

# ecoregions <- st_read("/PASANTIA/GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")

# smp_rect <- grid[subset(ecoregions, ECO_NAME == "Santa Marta páramo"),] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
# smp_rect$eco <- "SMP"

cmp_rect <- grid[subset(ecoregions, ECO_NAME == "Cordillera de Merida páramo"),] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
cmp_rect$eco <- "CMP"

ccp_rect <- grid[subset(ecoregions, ECO_NAME == "Cordillera Central páramo"),] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
ccp_rect$eco <- "CCP"

emm_rect <- grid[subset(ecoregions, ECO_NAME == "Ethiopian montane moorlands"),] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
emm_rect$eco <- "EMM"

eamm_rect <- grid[subset(ecoregions, ECO_NAME == "East African montane moorlands"),] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
eamm_rect$eco <- "EAMM"

crsg_rect <- grid[subset(ecoregions, ECO_NAME == "Central Range sub-alpine grasslands"),] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
crsg_rect$eco <- "CRSG"

nap <- subset(ecoregions, ecoregions$ECO_NAME == "Northern Andean páramo")
nap_rect <- grid[nap,] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
nap_rect$eco <- "NAP"

rects <- rbind(cmp_rect, ccp_rect, nap_rect, emm_rect, eamm_rect, crsg_rect)

# Rotate layers for png
rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)

eco.png.rot <- eco.png %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry")

box.png <- st_bbox(eco.png.rot)

world.rot <- world %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry")

ecoregions_dis.rot <- ecoregions %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry") %>%
  group_by(ECO_NAME, ECO_ID, G200_REGIO) %>%
  summarize()

crsg_rect.rot <- crsg_rect %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry")


#Create empty variables to store figures

maps <- list()
bp <- list()

#Produce figures iteratively
hspaces <- c("long_term", "extreme", "seasonality")

for (hspace in hspaces) {
  
  load(paste0("FIGURES/FINAL/", hspace, "_data.RData"))
  
  for (eco in 1:length(lpolygons.h)) {
    
    a <- rbind(lpolygons.h[[1]], lpolygons.h[[2]], lpolygons.h[[3]], lpolygons.h[[4]], lpolygons.h[[5]], lpolygons.h[[6]])
    a <- a[-eco,]
    b <- st_filter(lpoints.diff[[eco]], a)
    pixels.diff$intersect[which(pixels.diff$ID %in% b$ID)] <- -1
  }
  
  pixels.all <- rbind(pixels.diff, pixels.intersect)
  
  eco.grid.pols <- st_read("GIS/GRID/ECO_POL_INTERSECT.shp")
  eco.grid.pols <- eco.grid.pols[-c(1,5)]
  names(eco.grid.pols)[1:4] <- names(pixels.all)[1:4]
  
  eco.grid.pols$ecoregion[eco.grid.pols$ecoregion == "Santa Marta páramo"] <- "Northern Andean páramo"
  
  polis.all <- left_join(eco.grid.pols, st_drop_geometry(pixels.all), by = c("ecoregion", "G200_REGION", "polygon_id", "ID"))
  polis.all$intersect <- factor(polis.all$intersect, levels = c(1, -1, 0), ordered = TRUE)
  
  #Save polygons for summary figure
  saveRDS(polis.all, file = paste0("FIGURES/FINAL/", hspace, "_poligons.Rds"))
  
  # grid <- st_read("GIS/GRID/GRID.shp")
  
  #South america
  map.sam <- ggplot() +
    geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
    geom_sf(data = polis.all, aes(fill = intersect), color = NA, size = 0.01) +
    geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
    coord_sf(xlim = c(box.sam[1], box.sam[3]), ylim = c(box.sam[2], box.sam[4]), expand = TRUE) +
    annotate(geom = "text", x = -74, y = 2, label = "Colombia", size = 1.5, color = "grey22") +
    annotate(geom = "text", x = -72, y = -7.5, label = "Brazil", size = 1.5, color = "grey22") +
    annotate(geom = "text", x = -75, y = -5, label = "Peru", size = 1.5, color = "grey22") +
    annotate(geom = "text", x = -76.6, y = -0.75, label = "Ecuador", size = 1.5, color = "grey22") +
    annotate(geom = "text", x = -70.75, y = 7.4, label = "Venezuela", size = 1.5, color = "grey22") +
    annotate(geom = "text", x = -73.2, y = -3.5, label = "NAP", size = 2.5, color = "black", fontface = "bold") +
    annotate(geom = "text", x = -77.6, y = -4.5, label = "CCP", size = 2.5, color = "black", fontface = "bold") +
    annotate(geom = "text", x = -70.9, y = 10.1, label = "CMP", size = 2.5, color = "black", fontface = "bold") +
    scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00"), labels = c("Safe", "Analog", "Novel")) +
    scale_x_continuous(breaks = seq(-80, -70, 2),
                       labels = c(-80, -78, -76, -74, -72, -70)) +
    scale_y_continuous(breaks = seq(-10, 10, 5),
                       labels = c(-10, -5, 0, 5, 10)) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_legend(ticks.colour = "black",
                               frame.colour = "black",
                               title = NULL,
                               reverse = FALSE)) +
    annotation_scale(aes(location = "br",
                         style = "ticks"),
                     pad_x = unit(0.05, "cm"),
                     pad_y = unit(0.05, "cm"),
                     width_hint = 0.5,
                     height = unit(0.15, "cm"),
                     text_pad = unit(0.05, "cm"),
                     text_cex = 0.45) +
    annotation_north_arrow(location = "br",
                           pad_y = unit(0.5, "cm"),
                           height = unit(0.5, "cm"),
                           width = unit(0.5, "cm"),
                           style = north_arrow_fancy_orienteering(text_size = 4)) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          axis.text = element_text(size = 6),
          legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 10, face = "bold"))
  
  map.sam
  
  maps[[paste0(hspace, "_sam")]] <- map.sam
  
    # Africa
  map.afr <- ggplot() +
    geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
    geom_sf(data = polis.all, aes(fill = intersect), color = NA, size = 0.01) +
    geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
    coord_sf(xlim = c(box.afr[1], box.afr[3]), ylim = c(box.afr[2], box.afr[4]), expand = TRUE) +
    annotate(geom = "text", x = 39, y = 0, label = "Kenya", color = "grey22", size = 1.5) +
    annotate(geom = "text", x = 35, y = 8, label = "Ethiopia", color = "grey22", size = 1.5) +
    annotate(geom = "text", x = 35.5, y = -2.7, label = "Tanzania", color = "grey22", size = 1.5) +
    annotate(geom = "text", x = 36.5, y = 0.8, label = "EAMM", size = 2.5, color = "black", fontface = "bold") +
    annotate(geom = "text", x = 39.5, y = 5.9, label = "EMM", size = 2.5, color = "black", fontface = "bold") +
    scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00"), labels = c("Safe", "Analog", "Novel")) +
    scale_x_continuous(breaks = seq(35, 40, 2),
                       labels = c(35, 37, 39)) + 
    scale_y_continuous(breaks = seq(0, 10, 5),
                       labels = c(0, 5, 10)) + 
    labs(x = NULL, y = NULL) +
    annotation_scale(aes(location = "br",
                         style = "ticks"),
                     pad_x = unit(0.05, "cm"),
                     pad_y = unit(0.05, "cm"),
                     width_hint = 0.5,
                     height = unit(0.15, "cm"),
                     text_pad = unit(0.05, "cm"),
                     text_cex = 0.45) +
    annotation_north_arrow(location = "br",
                           pad_y = unit(0.5, "cm"),
                           height = unit(0.5, "cm"),
                           width = unit(0.5, "cm"),
                           style = north_arrow_fancy_orienteering(text_size = 4)) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          axis.text = element_text(size = 6),
          legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 10, face = "bold"))
  
  map.afr
  
  maps[[paste0(hspace, "_afr")]] <- map.afr
  
  # New Guinea
  
  polis.all.rot<- polis.all %>%
    mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
    st_drop_geometry() %>%
    rename(geometry = geom_rot) %>%
    st_set_geometry("geometry")
  
  map.png <- ggplot() +
    geom_sf(data = world.rot, color = gray(0.4), fill = gray(.93), size = 0.1) + 
    geom_sf(data = polis.all.rot, aes(fill = intersect), color = NA, size = 0.01) +
    geom_sf(data = crsg_rect.rot, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
    coord_sf(xlim = c(box.png[1], box.png[3]), ylim = c(box.png[2], box.png[4]), expand = TRUE) +
    scale_x_continuous(breaks = seq(-8, -4, 2),
                       labels = c(-8, -6, -4)) + 
    scale_y_continuous(breaks = seq(-148, -136, 4),
                       labels = c(148, 144, 140, 136)) +
    annotate(geom = "text", x = -6.5, y = -140, label = "Indonesia", size = 1.5, color = "grey22") +
    annotate(geom = "text", x = -7, y = -141.5, label = "Papua New Guinea", size = 1.5, color = "grey22") +
    annotate(geom = "text", x = -6, y = -135.8, label = "CRSG", size = 2.5, color = "black", fontface = "bold") +
    scale_fill_manual(values = c("#0072B2", "#E69F00", "#D55E00"), labels = c("Safe", "Analog", "Novel")) +
    labs(x = NULL, y = NULL) +
    guides(fill = guide_legend(ticks.colour = "black",
                               frame.colour = "black",
                               title = NULL,
                               reverse = FALSE)) +
    annotation_scale(aes(location = "br",
                         style = "ticks"),
                     plot_unit = "km",
                     pad_x = unit(0.05, "cm"),
                     pad_y = unit(0.05, "cm"),
                     width_hint = 0.5,
                     height = unit(0.15, "cm"),
                     text_pad = unit(0.05, "cm"),
                     text_cex = 0,
                     text_col = NA) +
    annotation_north_arrow(rotation = 270,
                           location = "br",
                           height = unit(0.5, "cm"),
                           width = unit(0.5, "cm"),
                           pad_y = unit(0.5, "cm"),
                           style = north_arrow_fancy_orienteering(
                             text_size = 4,
                             text_angle = 90
                           )) +
    annotate(geom = "text", x = -6.63, y = -148, label = "300 km", color = "black", size = 1.9) +
    theme_bw() +
    theme(panel.background = element_rect(fill = "white"),
          axis.text = element_text(size = 6),
          legend.position = "none",
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.title.y = element_text(size = 10, face = "bold"))
  
  map.png
  
  maps[[paste0(hspace, "_png")]] <- map.png
   
  # Barplots ----------------------------------------------------------------
  
  polis.all$area_km2 <- st_area(polis.all)
  units(polis.all$area_km2) <- "km2" #m2 to km2
  df.barplots <- st_drop_geometry(polis.all)
  df.barplots$condition <- ifelse(df.barplots$intersect == -1, "Analog",
                                  ifelse(df.barplots$intersect == 0, "Novel", "Safe"))
  
  df.barplots$area_km2 <- as.numeric(df.barplots$area_km2)
  
  df.barplots <- df.barplots %>%
    group_by(ecoregion, condition, intersect) %>%
    summarize(area_km2 = sum(area_km2))
  
  for (eco in unique(ecoregions$ECO_NAME)) {
    
    if (length(setdiff(c("Analog", "Novel", "Safe"), df.barplots$condition[df.barplots$ecoregion == eco])) > 0) {
      
      for (diff in setdiff(c("Analog", "Novel", "Safe"), df.barplots$condition[df.barplots$ecoregion == eco])) {
        
        df.barplots[nrow(df.barplots) + 1,] <- list(eco,
                                                    diff,
                                                    ifelse(diff == "Analog", "-1",
                                                           ifelse(diff == "Novel", "0", "1")),
                                                    0)
      }
    }
  }
  
  df.barplots.tae <- df.barplots %>%
    group_by(condition) %>%
    summarize(area_km2 = sum(area_km2))
  
  df.barplots.tae$ecoregion <- "Tropical Alpine Ecosystem"
  df.barplots <- rbind(df.barplots, df.barplots.tae)
  
  df.barplots <- df.barplots %>%
    group_by(ecoregion) %>%
    mutate(percent = (area_km2/sum(area_km2))*100)
  
  
  
  df.barplots$condition <- factor(df.barplots$condition, levels = c("Safe", "Analog", "Novel"), ordered = TRUE)
  df.barplots$ecoregion <- factor(df.barplots$ecoregion,
                                  levels = c("Northern Andean páramo",
                                             "Cordillera de Merida páramo",
                                             "Cordillera Central páramo",
                                             "East African montane moorlands",
                                             "Ethiopian montane moorlands",
                                             "Central Range sub-alpine grasslands",
                                             "Tropical Alpine Ecosystem"))
  
  #Save df for summary figure
  saveRDS(df.barplots, file = paste0("FIGURES/FINAL/", hspace, "_df.Rds"))
  
  barplot <- ggplot(df.barplots, aes(fill = condition, x = ecoregion, y = percent)) +
    geom_bar(position = "dodge", stat = "identity", width = 0.8, colour = "black") +
    xlab(NULL) +
    ylab("Percentage of area (%)") +
    scale_y_continuous(breaks = seq(0,100, 25),limits = c(0, 100),
                       expand = expansion(mult = c(0, 0.009))) +
    scale_x_discrete(labels = c("NAP","CMP", "CCP", "EAMM", "EMM", "CRSG", "HETEs")) +
    scale_fill_manual(values = c("#0072B2BF", "#E69F00BF", "#D55E00BF")) +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(colour = "Black"))
  
  barplot
  
  bp[[hspace]] <- barplot
  
}

# Join figures

map <- (maps$long_term_sam + ylab("LT") | maps$long_term_afr + guides(fill = "none") | maps$long_term_png + guides(fill = "none")) /
  (maps$extreme_sam + ylab("EX") + guides(fill = "none") | maps$extreme_afr + guides(fill = "none") | maps$extreme_png + guides(fill = "none")) /
  (maps$seasonality_sam + ylab("SE") + guides(fill = "none") | maps$seasonality_afr + guides(fill = "none") | maps$seasonality_png + guides(fill = "none")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "", "", "B", "", "", "C", "", ""))) &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", hjust = 1, vjust = 3),
        axis.title.y = element_text(size = 10, face = "bold"))

map

ggsave(plot = map, "FIGURES/FINAL/FIGURE_S4_NEW.pdf", device = "pdf", dpi = 1000, width = 18, height = 27, units = "cm")
ggsave(plot = map, "FIGURES/FINAL/FIGURE_S4_NEW.svg", device = "svg", dpi = 1000, width = 18, height = 27, units = "cm")
ggsave(plot = map, "FIGURES/FINAL/FIGURE_S4_NEW.png", device = "png", dpi = 1000, width = 18, height = 27, units = "cm")


plot <- ((bp$long_term + labs(title = "LT")) / (bp$extreme + labs(title = "EX") + guides(fill = "none")) / (bp$seasonality + labs(title = "SE") + guides(fill = "none"))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "B", "C"))) &
  theme(legend.position = "bottom",
        plot.tag = element_text(face = "bold", hjust = 1, vjust = 3, size = 10),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.title = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text = element_text(size = 8))
  

plot

ggsave(plot = plot, "FIGURES/FINAL/FIGURE_S5_NEW.pdf", device = "pdf", dpi = 1000, width = 12, height = 15, units = "cm")
ggsave(plot = plot, "FIGURES/FINAL/FIGURE_S5_NEW.svg", device = "svg", dpi = 1000, width = 12, height = 15, units = "cm")
ggsave(plot = plot, "FIGURES/FINAL/FIGURE_S5_NEW.png", device = "png", dpi = 1000, width = 12, height = 15, units = "cm")


# Safe hydroclimatic space changes summary --------------------------------

col <- rev(c("#238b45", "yellow", "#FFCC00", "#FF9900", "#FF6600", "#FF3300", "#CC0000", "darkred"))
#Prepare TAE polygons
lt.pol <- readRDS(file = "FIGURES/FINAL/long_term_poligons.Rds")
ex.pol <- readRDS(file = "FIGURES/FINAL/extreme_poligons.Rds")
se.pol <- readRDS(file = "FIGURES/FINAL/seasonality_poligons.Rds")
pol <- lt.pol

pol$intersect.lt <- ifelse(lt.pol$intersect == "0" | lt.pol$intersect == "-1", "1", "0")
pol$intersect.ex <- ifelse(ex.pol$intersect == "0" | ex.pol$intersect == "-1", "1", "0")
pol$intersect.se <- ifelse(se.pol$intersect == "0" | se.pol$intersect == "-1", "1", "0")

my_cols <- c("intersect.lt", "intersect.ex", "intersect.se")

pol$intersect.all <- do.call(paste, c(st_drop_geometry(pol[my_cols]), sep = ""))

pol$intersect.all <- factor(pol$intersect.all,
                                    levels = c("111", "110", "101", "011", "100", "010", "001", "000"), ordered = TRUE)


#Prepare spatial data
ecoregions <- st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")
ecoregions[which(ecoregions$ECO_NAME == "Santa Marta páramo"),][4] <- "Northern Andean páramo"
ecoregions_dis <- ecoregions %>% group_by(ECO_NAME, G200_REGIO) %>% summarize()
ecoregions_dis$ID <- seq(1:nrow(ecoregions_dis))

eco.sam <- subset(ecoregions, ecoregions$ECO_NAME == "Cordillera Central páramo" |
                    ecoregions$ECO_NAME == "Northern Andean páramo" |
                    ecoregions$ECO_NAME == "Cordillera de Merida páramo")

box.sam <- st_bbox(eco.sam)

eco.afr <- subset(ecoregions, ecoregions$ECO_NAME == "East African montane moorlands" |
                    ecoregions$ECO_NAME == "Ethiopian montane moorlands")

box.afr <- st_bbox(eco.afr)

eco.png <- subset(ecoregions, ecoregions$ECO_NAME == "Central Range sub-alpine grasslands")

world <- ne_countries(scale = "medium", returnclass = "sf")

# Rotate layers for png
rot <- function(a) matrix(c(cos(a), sin(a), -sin(a), cos(a)), 2, 2)

eco.png.rot <- eco.png %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry")

box.png <- st_bbox(eco.png.rot)

world.rot <- world %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry")

ecoregions_dis.rot <- ecoregions %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry") %>%
  group_by(ECO_NAME, ECO_ID, G200_REGIO) %>%
  summarize()

#South america
map.sam <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) +
  geom_sf(data = pol, aes(fill = intersect.all), color = NA) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.sam[1], box.sam[3]), ylim = c(box.sam[2], box.sam[4]), expand = TRUE) +
  annotate(geom = "text", x = -74, y = 2, label = "Colombia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -72, y = -7.5, label = "Brazil", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -75, y = -5, label = "Peru", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -76.6, y = -0.75, label = "Ecuador", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -70.75, y = 7.4, label = "Venezuela", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -73.2, y = -3.5, label = "NAP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -77.6, y = -4.5, label = "CCP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -70.9, y = 10.1, label = "CMP", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_manual(values = col, labels = c("LT + EX + SE", "LT + EX", "LT + SE", "EX + SE", "LT", "EX", "SE", "None")) +
  scale_x_continuous(breaks = seq(-80, -70, 2),
                     labels = c(-80, -78, -76, -74, -72, -70)) +
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = c(-10, -5, 0, 5, 10)) +
  labs(x = NULL, y = NULL, title = "Neotropic") +
  guides(fill = guide_legend(ticks.colour = "black",
                             frame.colour = "black",
                             title = "Breached spaces",
                             reverse = TRUE,
                             nrow = 1,
                             title.position = "top",
                             title.theme = element_text(hjust = 0.5, size = 10))) +
  annotation_scale(aes(location = "br",
                       style = "ticks"),
                   pad_x = unit(0.05, "cm"),
                   pad_y = unit(0.05, "cm"),
                   width_hint = 0.5,
                   height = unit(0.15, "cm"),
                   text_pad = unit(0.05, "cm"),
                   text_cex = 0.45) +
  annotation_north_arrow(location = "br",
                         pad_y = unit(0.5, "cm"),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 4)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

map.sam

# Africa
map.afr <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_sf(data = pol, aes(fill = intersect.all), color = NA, size = 0.01) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.afr[1], box.afr[3]), ylim = c(box.afr[2], box.afr[4]), expand = TRUE) +
  annotate(geom = "text", x = 39, y = 0, label = "Kenya", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35, y = 8, label = "Ethiopia", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35.5, y = -2.7, label = "Tanzania", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 36.5, y = 0.8, label = "EAMM", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = 39.5, y = 5.9, label = "EMM", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_manual(values = col, labels = c("LT + EX + SE", "LT + EX", "LT + SE", "EX + SE", "LT", "EX", "SE", "None")) +
  scale_x_continuous(breaks = seq(35, 40, 2),
                     labels = c(35, 37, 39)) + 
  scale_y_continuous(breaks = seq(0, 10, 5),
                     labels = c(0, 5, 10)) + 
  labs(x = NULL, y = NULL, title = "Afrotropic") +
  guides(fill = guide_legend(ticks.colour = "black",
                             frame.colour = "black",
                             title = "Breached spaces",
                             reverse = TRUE,
                             nrow = 1,
                             title.position = "top",
                             title.theme = element_text(hjust = 0.5, size = 10))) +
  annotation_scale(aes(location = "br",
                       style = "ticks"),
                   pad_x = unit(0.05, "cm"),
                   pad_y = unit(0.05, "cm"),
                   width_hint = 0.5,
                   height = unit(0.15, "cm"),
                   text_pad = unit(0.05, "cm"),
                   text_cex = 0.45) +
  annotation_north_arrow(location = "br",
                         pad_y = unit(0.5, "cm"),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 4)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

map.afr

# New Guinea

pol.rot<- pol %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry")

map.png <- ggplot() +
  geom_sf(data = world.rot, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_sf(data = pol.rot, aes(fill = intersect.all), color = NA, size = 0.01) +
  geom_sf(data = crsg_rect.rot, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.png[1], box.png[3]), ylim = c(box.png[2], box.png[4]), expand = TRUE) +
  scale_x_continuous(breaks = seq(-8, -4, 2),
                     labels = c(-8, -6, -4)) + 
  scale_y_continuous(breaks = seq(-148, -136, 4),
                     labels = c(148, 144, 140, 136)) +
  annotate(geom = "text", x = -6.5, y = -140, label = "Indonesia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -7, y = -141.5, label = "Papua New Guinea", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -6, y = -135.8, label = "CRSG", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_manual(values = col, labels = c("LT + EX + SE", "LT + EX", "LT + SE", "EX + SE", "LT", "EX", "SE", "None")) +
  labs(x = NULL, y = NULL, title = "Australasia") +
  guides(fill = guide_legend(ticks.colour = "black",
                             frame.colour = "black",
                             title = "Breached spaces",
                             reverse = TRUE,
                             nrow = 1, 
                             title.position = "top",
                             title.theme = element_text(hjust = 0.5, size = 10))) +
  annotation_scale(aes(location = "br",
                       style = "ticks"),
                   plot_unit = "km",
                   pad_x = unit(0.05, "cm"),
                   pad_y = unit(0.05, "cm"),
                   width_hint = 0.5,
                   height = unit(0.15, "cm"),
                   text_pad = unit(0.05, "cm"),
                   text_cex = 0,
                   text_col = NA) +
  annotation_north_arrow(rotation = 270,
                         location = "br",
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         pad_y = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(
                           text_size = 4,
                           text_angle = 90
                         )) +
  annotate(geom = "text", x = -6.63, y = -148, label = "300 km", color = "black", size = 1.9) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        legend.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5))

map.png

# Summary Barplots ----------------------------------------------------------------

pol$area_km2 <- st_area(pol)
units(pol$area_km2) <- "km2" #m2 to km2
df.barplots <- st_drop_geometry(pol)

df.barplots$area_km2 <- as.numeric(df.barplots$area_km2)

df.barplots <- df.barplots %>%
  group_by(ecoregion, intersect.all) %>%
  summarize(area_km2 = sum(area_km2))

for (eco in unique(ecoregions$ECO_NAME)) {
  
  if (length(setdiff(c("111", "101", "001", "000", "100", "110", "010", "011"), df.barplots$intersect.all[df.barplots$ecoregion == eco])) > 0) {
    
    for (diff in setdiff(c("111", "101", "001", "000", "100", "110", "010", "011"), df.barplots$intersect.all[df.barplots$ecoregion == eco])) {
      
      df.barplots[nrow(df.barplots) + 1,] <- list(eco, diff, 0)
    }
  }
}

df.barplots.tae <- df.barplots %>%
  group_by(intersect.all) %>%
  summarize(area_km2 = sum(area_km2))

df.barplots.tae$ecoregion <- "Tropical Alpine Ecosystem"
df.barplots <- rbind(df.barplots, df.barplots.tae)

df.barplots <- df.barplots %>%
  group_by(ecoregion) %>%
  mutate(percent = (area_km2/sum(area_km2))*100)



df.barplots$intersect.all <- factor(df.barplots$intersect.all,
                                    levels = c("111", "110", "101", "011", "100", "010", "001", "000"), ordered = TRUE)

df.barplots$ecoregion <- factor(df.barplots$ecoregion,
                                levels = c("Northern Andean páramo",
                                           "Cordillera de Merida páramo",
                                           "Cordillera Central páramo",
                                           "East African montane moorlands",
                                           "Ethiopian montane moorlands",
                                           "Central Range sub-alpine grasslands",
                                           "Tropical Alpine Ecosystem"))

df.barplots$category <- ifelse(df.barplots$intersect.all %in% c("000"), "Low",
          ifelse(df.barplots$intersect.all %in% c("001", "010", "100"), "Medium",
                 ifelse(df.barplots$intersect.all %in% c("011", "101", "110"), "High",
                        ifelse(df.barplots$intersect.all %in% c("111"), "Critical", NA))))

df.barplots$category <- factor(df.barplots$category,
                      levels = c("Low", "Medium", "High", "Critical"), ordered = TRUE)

df.barplots$biorealm <- ifelse(df.barplots$ecoregion %in% c("Santa Marta páramo", "Northern Andean páramo",
                                                 "Cordillera de Merida páramo", "Cordillera Central páramo"), "Neotropic",
                          ifelse(df.barplots$ecoregion %in% c("East African montane moorlands", "Ethiopian montane moorlands"), "Afrotropic",
                                 ifelse(df.barplots$ecoregion == "Central Range sub-alpine grasslands", "Australasia", "")))

df.barplots <- df.barplots %>%
  mutate(biorealm = fct_relevel(biorealm, "Neotropic", "Afrotropic", "Australasia", ""))

df.barplots$eco_acr <- ifelse(df.barplots$ecoregion == "Santa Marta páramo", "SMP",
                         ifelse(df.barplots$ecoregion == "Northern Andean páramo", "NAP",
                                ifelse(df.barplots$ecoregion == "Cordillera de Merida páramo", "CMP",
                                       ifelse(df.barplots$ecoregion == "Cordillera Central páramo", "CCP",
                                              ifelse(df.barplots$ecoregion == "East African montane moorlands", "EAMM",
                                                     ifelse(df.barplots$ecoregion == "Ethiopian montane moorlands", "EMM",
                                                            ifelse(df.barplots$ecoregion == "Central Range sub-alpine grasslands", "CRSG",
                                                                   "HETEs")))))))

df.barplots <- df.barplots %>%
  mutate(eco_acr = fct_relevel(eco_acr, "SMP", "NAP", "CMP", "CCP", "EAMM", "EMM", "CRSG", "HETEs"))

#With hatch
barplot <- ggplot(df.barplots, aes(fill = intersect.all, x = interaction(eco_acr, biorealm, sep = "!"), y = percent, pattern = category)) +
  geom_bar_pattern(position = "stack", stat = "identity", width = 0.7,
                   color = NA, 
                   pattern_fill = "black",
                   pattern_density = 0.1,
                   pattern_spacing = 0.025,
                   pattern_angle = 45,
                   pattern_key_scale_factor = 0.6,
                   pattern_size = 0.1) +
  geom_vline(xintercept = 6.5, linetype = "dashed", size = 0.4) +
  labs(x = NULL, y = "Percentage of area (%)", pattern = "Vulnerability", fill = "Breached spaces") +
  scale_y_continuous(n.breaks = 6, expand = expansion(mult = c(0, 0.009))) +
  scale_x_discrete(guide = guide_axis_nested(delim = "!", extend = 0.75)) +
  scale_pattern_manual(values = c("none", "stripe", "circle", "crosshatch")) +
  scale_fill_manual(values = col) +
  guides(fill = guide_legend(override.aes = list(pattern = "none"), nrow = 1, title.position = "top", hjust = 0.5),
         pattern = guide_legend(override.aes = list(fill = "white", color = "black"), title.position = "top", hjust = 0.5)) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(hjust = 0.5, size = 10),
        legend.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        axis.text = element_text(colour = "black", size = 8),
        axis.ticks = element_line(colour = "black"),
        ggh4x.axis.nesttext.x = element_text(size = 8),
        ggh4x.axis.nestline.x = element_line(colour = "black"))

barplot

# Join figures

layout <- "
ABC
DDD"

map <- ((map.sam | map.afr  | map.png) + barplot + guides(fill = "none")) +
  plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = list(c("A", "", "", "B"))) &
  theme(legend.position = "bottom", legend.box = "vertical", plot.tag = element_text(face = "bold", hjust = 1, vjust = 1))

map

ggsave(plot = map, "FIGURES/FINAL/FIGURE4_NEW.pdf", device = "pdf", dpi = 1000, width = 160, height = 200, units = "mm")
ggsave(plot = map, "FIGURES/FINAL/FIGURE4_NEW.svg", device = "svg", dpi = 1000, width = 160, height = 200, units = "mm")
ggsave(plot = map, "FIGURES/FINAL/FIGURE4_NEW.png", device = "png", dpi = 1000, width = 160, height = 200, units = "mm")
