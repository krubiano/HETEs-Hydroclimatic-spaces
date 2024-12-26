##CREATE FIGURE 1##

library(stringr)
library(raster)
library(terra)
library(tidyterra)
library(sf)
library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(rnaturalearth)
library(ggspatial)
library(ggmap)
library(Cairo)


input <- "4_ANNUAL_STACKS_SUMMARY"

# Aggregates the 30-year rasters by model to produce an ensemble ----------

if (!dir.exists(file.path(input, "ENSEMBLE", "biovar1", "diff"))) {
  dir.create(file.path(input, "ENSEMBLE", "biovar1", "diff"), recursive = TRUE)}

if (!dir.exists(file.path(input, "ENSEMBLE", "biovar12", "diff"))) {
  dir.create(file.path(input, "ENSEMBLE", "biovar12", "diff"), recursive = TRUE)}

#Creates and filters  a list with file paths of tif files in the input directory
paths <- list.files(file.path(input), full.names = TRUE, recursive = TRUE, pattern = "\\.tif$")
paths <- paths[str_detect(paths, paste0("\\b", "diff", "\\b"))]
paths <- paths[str_detect(paths, paste0("\\b", "ENSEMBLE", "\\b"), negate = TRUE)]

#Stack rasters for tas and aggregate them to create the ensemble
paths.tas <- paths[str_detect(paths, paste0("\\b", "biovar1", "\\b"))]

r.tas <- stack(paths.tas)
r.tas.ens <- stackApply(r.tas,
                        indices = rep(1, 6),
                        fun = mean,
                        na.rm = TRUE,
                        filename = file.path(input, "ENSEMBLE", "biovar1", "diff", "biovar1_ens_diff.tif"),
                        format = "GTiff",
                        datatype = "FLT4S",
                        overwrite = TRUE)

r.tas.ens.err <- stackApply(r.tas,
                        indices = rep(1, 6),
                        fun = function(x, ...) sd(x)/sqrt(6),
                        filename = file.path(input, "ENSEMBLE", "biovar1", "diff", "biovar1_ens_diff_stderr.tif"),
                        format = "GTiff",
                        datatype = "FLT4S",
                        overwrite = TRUE)

#Stack rasters for pr and aggregate them to create the ensemble
paths.pr <- paths[str_detect(paths, paste0("\\b", "biovar12", "\\b"))]

r.pr <- stack(paths.pr)
r.pr.ens <- stackApply(r.pr,
                        indices = rep(1, 6),
                        fun = mean,
                        na.rm = TRUE,
                        filename = file.path(input, "ENSEMBLE", "biovar12", "diff", "biovar12_ens_diff.tif"),
                        format = "GTiff",
                        datatype = "FLT4S",
                        overwrite = TRUE)

r.pr.ens.err <- stackApply(r.pr,
                            indices = rep(1, 6),
                            fun = function(x, ...) sd(x)/sqrt(6),
                            filename = file.path(input, "ENSEMBLE", "biovar12", "diff", "biovar12_ens_diff_stderr.tif"),
                            format = "GTiff",
                            datatype = "FLT4S",
                            overwrite = TRUE)

### Load ecoregions shapefile and subset by location ------------------------

ecoregions <- st_read("/PASANTIA/GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")

eco.sam <- subset(ecoregions, ecoregions$ECO_NAME == "Cordillera Central páramo" |
                    ecoregions$ECO_NAME == "Cordillera de Merida páramo" |
                    ecoregions$ECO_NAME == "Northern Andean páramo" |
                    ecoregions$ECO_NAME == "Santa Marta páramo")

box.sam <- st_bbox(eco.sam)

eco.afr <- subset(ecoregions, ecoregions$ECO_NAME == "East African montane moorlands" |
                    ecoregions$ECO_NAME == "Ethiopian montane moorlands")

box.afr <- st_bbox(eco.afr)

eco.png <- subset(ecoregions, ecoregions$ECO_NAME == "Central Range sub-alpine grasslands")

world <- ne_countries(scale = "medium", returnclass = "sf")


# Prepare polygons to highlight ecoregions --------------------------------

grid <- st_read("GIS/GRID/GRID.shp")

smp_rect <- grid[subset(ecoregions, ECO_NAME == "Santa Marta páramo"),] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
smp_rect$eco <- "SMP"
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

# gmba <- st_read("GIS/GMBA_Inventory_v2.0_standard/GMBA_Inventory_v2.0_standard.shp")
# cordor <- subset(gmba, gmba$GMBA_V2_ID == 11265)

nap <- subset(ecoregions, ecoregions$ECO_NAME == "Northern Andean páramo")

nap_rect <- grid[nap,] %>% summarize %>% st_bbox() %>% st_as_sfc() %>% st_sf()
nap_rect$eco <- "NAP"

rects <- rbind(smp_rect, cmp_rect, ccp_rect, nap_rect, emm_rect, eamm_rect, crsg_rect)

### Rotate Papua New Guinea shapefile  --------------------------------------

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


crsg_rect.rot <- crsg_rect %>%
  mutate(geom_rot = st_geometry(.)*rot(pi/2)) %>%
  st_drop_geometry() %>%
  rename(geometry = geom_rot) %>%
  st_set_geometry("geometry")

# Create inset ------------------------------------------------------------

inset <- ggplot() +
  geom_sf(data = world, fill = gray(0.7), color = NA) +
  layer_spatial(data = box.sam, fill = NA, color = "red", linewidth = 0.25) +
  layer_spatial(data = box.afr, fill = NA, color = "red", linewidth = 0.25) +
  layer_spatial(data = st_bbox(eco.png), fill = NA, color = "red", size = 0.25) +
  theme_classic() +
  theme(plot.margin = margin(0),
        axis.title = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black"),
        plot.background = element_rect(fill = "aliceblue"))

inset

### Difference maps ----------------------------------------------

## Precipitation -----------------------------------------------------------

pr.diff <- rast(file.path(input, "ENSEMBLE", "biovar12", "diff", "biovar12_ens_diff.tif"))

limit <- max(abs(values(pr.diff, na.rm = TRUE))) * c(-1, 1)

pr.diff.se <- rast(file.path(input, "ENSEMBLE", "biovar12", "diff", "biovar12_ens_diff_stderr.tif"))

# South America -----------------------------------------------------------

#Mean
pr.diff.sam.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = pr.diff, maxcell = Inf) +
  geom_sf(data = st_simplify(ecoregions, preserveTopology = TRUE, dTolerance = 1000), color = gray(.5), fill = NA, size = 0.01) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.sam[1], box.sam[3]), ylim = c(box.sam[2], box.sam[4]), expand = TRUE) +
  annotate(geom = "text", x = -74, y = 2, label = "Colombia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -72, y = -7.5, label = "Brazil", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -75, y = -5, label = "Peru", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -76.6, y = -0.75, label = "Ecuador", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -70.75, y = 7.4, label = "Venezuela", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -73.2, y = -3.5, label = "NAP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -77.6, y = -4.5, label = "CCP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -73.6, y = 11.5, label = "SMP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -70.9, y = 10.1, label = "CMP", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_gradientn(colours = c("darkred", "red", "orange", "yellow", "white", "skyblue", "royalblue", "mediumblue", "darkblue"), na.value = NA,
                       limits = limit) +
  scale_x_continuous(breaks = seq(-80, -70, 2),
                     labels = c(-80, -78, -76, -74, -72, -70)) + 
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = c(-10, -5, 0, 5, 10)) + 
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = paste0("\u0394", "P", " (mm)"),
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

pr.diff.sam.fig

# sam.inset <- ggplot() +
#   geom_sf(data = world, fill = gray(0.7), color = NA) +
#   # geom_hline(yintercept = 0, color = gray(0.6), linewidth = 0.2, linetype = "solid") +
#   # geom_hline(yintercept = 23.27, color = gray(0.6), linewidth = 0.2, linetype = "dashed") +
#   # geom_hline(yintercept = -23.27, color = gray(0.6), linewidth = 0.2, linetype = "dashed") +
#   layer_spatial(data = box.sam, fill = NA, color = "red", linewidth = 0.25) +
#   theme_classic() +
#   theme(plot.margin = margin(0),
#         axis.title = element_blank(),
#         panel.border = element_rect(fill = NA, colour = "black"),
#         plot.background = element_rect(fill = "aliceblue"))
# 
# sam.inset
# 
# sam.x.range <- ggplot_build(pr.diff.sam.plot)$layout$panel_params[[1]]$x_range
# sam.y.range <- ggplot_build(pr.diff.sam.plot)$layout$panel_params[[1]]$y_range
# 
# pr.diff.sam.fig <- pr.diff.sam.plot + annotation_custom(ggplotGrob(sam.inset), xmin = sam.x.range[1], xmax = sam.x.range[1]+4, ymax = sam.y.range[2], ymin = sam.y.range[2]-2)
# 
# pr.diff.sam.fig

#Standard error
pr.diff.se.sam.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = pr.diff.se, maxcell = Inf) +
  geom_sf(data = st_simplify(ecoregions, preserveTopology = TRUE, dTolerance = 1000), color = gray(.5), fill = NA, size = 0.01) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.sam[1], box.sam[3]), ylim = c(box.sam[2], box.sam[4]), expand = TRUE) +
  annotate(geom = "text", x = -74, y = 2, label = "Colombia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -72, y = -7.5, label = "Brazil", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -75, y = -5, label = "Peru", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -76.6, y = -0.75, label = "Ecuador", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -70.75, y = 7.4, label = "Venezuela", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -73.2, y = -3.5, label = "NAP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -77.6, y = -4.5, label = "CCP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -73.6, y = 11.5, label = "SMP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -70.9, y = 10.1, label = "CMP", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_gradientn(colours = c("white", "skyblue", "royalblue", "mediumblue", "darkblue"), na.value = NA,
                       limits = c(0, 160)) +
  scale_x_continuous(breaks = seq(-80, -70, 2),
                     labels = c(-80, -78, -76, -74, -72, -70)) + 
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = c(-10, -5, 0, 5, 10)) + 
  labs(x = NULL, y = "Precipitation") +
  annotation_scale(aes(location = "br",
                       style = "ticks"),
                   pad_x = unit(0.05, "cm"),
                   pad_y = unit(0.05, "cm"),
                   width_hint = 0.5,
                   height = unit(0.15, "cm"),
                   text_pad = unit(0.05, "cm"),
                   text_cex = 0.45) +
  annotation_north_arrow(location = "br",
                         # pad_x = unit(0, "cm"),
                         pad_y = unit(0.5, "cm"),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 4)) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = "mm",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

pr.diff.se.sam.fig


# pr.diff.se.sam.fig <- pr.diff.se.sam.plot + annotation_custom(ggplotGrob(sam.inset), xmin = sam.x.range[1], xmax = sam.x.range[1]+4, ymax = sam.y.range[2], ymin = sam.y.range[2]-2)
# pr.diff.se.sam.fig 

# Africa ------------------------------------------------------------------

#mean
pr.diff.afr.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = pr.diff, maxcell = Inf) +
  geom_sf(data = ecoregions, color = gray(.5), fill = NA, size = 0.15) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.afr[1], box.afr[3]), ylim = c(box.afr[2], box.afr[4]), expand = TRUE) +
  annotate(geom = "text", x = 39, y = 0, label = "Kenya", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35, y = 8, label = "Ethiopia", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35.5, y = -2.7, label = "Tanzania", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 36.5, y = 0.8, label = "EAMM", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = 39.5, y = 5.9, label = "EMM", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_gradientn(colours = c("darkred", "red", "orange", "yellow", "white", "skyblue", "royalblue", "mediumblue", "darkblue"),na.value = NA,
                       limits = limit) +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = paste0("\u0394", "P", " (mm)"),
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

pr.diff.afr.fig

# afr.inset <- ggplot() +
#   geom_sf(data = world, fill = gray(0.7), color = NA) +
#   # geom_hline(yintercept = 0, color = gray(0.6), size = 0.2, linetype = "solid") +
#   # geom_hline(yintercept = 23.27, color = gray(0.6), size = 0.2, linetype = "dashed") +
#   # geom_hline(yintercept = -23.27, color = gray(0.6), size = 0.2, linetype = "dashed") +
#   layer_spatial(data = box.afr, fill = NA, color = "red", size = 0.25) +
#   theme_classic() +
#   theme(plot.margin = margin(0),
#         axis.title = element_blank(),
#         panel.border = element_rect(fill = NA, colour = "black"),
#         plot.background = element_rect(fill = "aliceblue"))
# 
# afr.inset
# 
# afr.x.range <- ggplot_build(pr.diff.afr.plot)$layout$panel_params[[1]]$x_range
# afr.y.range <- ggplot_build(pr.diff.afr.plot)$layout$panel_params[[1]]$y_range
# 
# pr.diff.afr.fig <- pr.diff.afr.plot + annotation_custom(ggplotGrob(afr.inset), xmin = afr.x.range[1], xmax = afr.x.range[1]+3.31, ymax = afr.y.range[2], ymin = afr.y.range[2]-1.59)
# 
# pr.diff.afr.fig

#Standard error
pr.diff.se.afr.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = pr.diff.se, maxcell = Inf) +
  geom_sf(data = ecoregions, color = gray(.5), fill = NA, size = 0.01) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.afr[1], box.afr[3]), ylim = c(box.afr[2], box.afr[4]), expand = TRUE) +
  annotate(geom = "text", x = 39, y = 0, label = "Kenya", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35, y = 8, label = "Ethiopia", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35.5, y = -2.7, label = "Tanzania", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 36.5, y = 0.8, label = "EAMM", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = 39.5, y = 5.9, label = "EMM", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_gradientn(colours = c("white", "skyblue", "royalblue", "mediumblue", "darkblue"), na.value = NA,
                       limit = c(0, 160)) +
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
                         # pad_x = unit(0, "cm"),
                         pad_y = unit(0.5, "cm"),
                         height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"),
                         style = north_arrow_fancy_orienteering(text_size = 4)) +
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = "mm",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

pr.diff.se.afr.fig

# pr.diff.se.afr.fig <- pr.diff.se.afr.plot + annotation_custom(ggplotGrob(afr.inset), xmin = afr.x.range[1], xmax = afr.x.range[1]+3.31, ymax = afr.y.range[2], ymin = afr.y.range[2]-1.59)
# pr.diff.se.afr.fig 


# Papua New Guinea --------------------------------------------------------

pr.diff.r <- t(flip(stack(pr.diff), 2))
pr.diff.r.df <- as.data.frame(pr.diff.r, xy = TRUE) %>% drop_na()
names(pr.diff.r.df) <- c("x", "y", "X1")

pr.diff.se.r <- t(flip(stack(pr.diff.se), 2))
pr.diff.se.r.df <- as.data.frame(pr.diff.se.r, xy = TRUE) %>% drop_na()
names(pr.diff.se.r.df) <- c("x", "y", "X1")

#mean
pr.diff.png.fig <- ggplot() +
  geom_sf(data = world.rot, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_raster(aes(x = x, y = y, fill = X1), data = pr.diff.r.df) +
  geom_sf(data = crsg_rect.rot, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  geom_sf(data = eco.png.rot, color = gray(.5), fill = NA, size = 0.01) +
  scale_x_continuous(breaks = seq(-8, -4, 2),
                     labels = c(-8, -6, -4)) + 
  scale_y_continuous(breaks = seq(-148, -136, 4),
                     labels = c(148, 144, 140, 136)) +
  annotate(geom = "text", x = -6.5, y = -140, label = "Indonesia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -7, y = -141.5, label = "Papua New Guinea", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -6.85, y = -148, label = "300 km", color = "black", size = 1.9) +
  annotate(geom = "text", x = -6, y = -135.8, label = "CRSG", size = 2.5, color = "black", fontface = "bold") +
  coord_sf(xlim = c(box.png[1], box.png[3]), ylim = c(box.png[2], box.png[4]), expand = TRUE) +
  scale_fill_gradientn(colours = c("darkred", "red", "orange", "yellow", "white", "skyblue", "royalblue", "mediumblue", "darkblue"),
                       limits = limit) +
  labs(x = NULL, y = NULL) +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = paste0("\u0394", "P", " (mm)"),
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

pr.diff.png.fig

# png.inset <- ggplot() +
#   geom_sf(data = world, fill = gray(0.7), color = NA) +
#   # geom_hline(yintercept = 0, color = gray(0.6), size = 0.2, linetype = "solid") +
#   # geom_hline(yintercept = 23.27, color = gray(0.6), size = 0.2, linetype = "dashed") +
#   # geom_hline(yintercept = -23.27, color = gray(0.6), size = 0.2, linetype = "dashed") +
#   layer_spatial(data = st_bbox(eco.png), fill = NA, color = "red", size = 0.25) +
#   theme_classic() +
#   theme(plot.margin = margin(0),
#         axis.title = element_blank(),
#         panel.border = element_rect(fill = NA, colour = "black"),
#         plot.background = element_rect(fill = "aliceblue"))
# 
# png.inset
# 
# png.x.range <- ggplot_build(pr.diff.png.plot)$layout$panel_params[[1]]$x_range
# png.y.range <- ggplot_build(pr.diff.png.plot)$layout$panel_params[[1]]$y_range
# 
# pr.diff.png.fig <- pr.diff.png.plot + annotation_custom(ggplotGrob(png.inset), xmin = png.x.range[1], xmax = png.x.range[1]+2.3, ymax = png.y.range[2], ymin = -136.2)
# 
# pr.diff.png.fig

#Standard error
pr.diff.se.png.fig <- ggplot() +
  geom_sf(data = world.rot, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_raster(aes(x = x, y = y, fill = X1), data = pr.diff.se.r.df) +
  geom_sf(data = crsg_rect.rot, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  geom_sf(data = eco.png.rot, color = gray(.5), fill = NA, size = 0.01) +
  scale_x_continuous(breaks = seq(-8, -4, 2),
                     labels = c(-8, -6, -4)) + 
  scale_y_continuous(breaks = seq(-148, -136, 4),
                     labels = c(148, 144, 140, 136)) +
  annotate(geom = "text", x = -6.5, y = -140, label = "Indonesia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -7, y = -141.5, label = "Papua New Guinea", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -6.85, y = -148, label = "300 km", color = "black", size = 1.9) +
  annotate(geom = "text", x = -6, y = -135.8, label = "CRSG", size = 2.5, color = "black", fontface = "bold") +
  coord_sf(xlim = c(box.png[1], box.png[3]), ylim = c(box.png[2], box.png[4]), expand = TRUE) +
  scale_fill_gradientn(colours = c("white", "skyblue", "royalblue", "mediumblue", "darkblue"),
                       limit = c(0, 160)) +
  labs(x = NULL, y = NULL) +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = "mm",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

pr.diff.se.png.fig

# pr.diff.se.png.fig <- pr.diff.se.png.plot + annotation_custom(ggplotGrob(png.inset), xmin = png.x.range[1], xmax = png.x.range[1]+2.3, ymax = png.y.range[2], ymin = -136.2)
# pr.diff.se.png.fig

## Temperature -----------------------------------------------------------

# #mean
tas.diff <- rast(file.path(input, "ENSEMBLE", "biovar1", "diff", "biovar1_ens_diff.tif"))

#standard error
tas.diff.se <- rast(file.path(input, "ENSEMBLE", "biovar1", "diff", "biovar1_ens_diff_stderr.tif"))

# South America -----------------------------------------------------------

#mean
tas.diff.sam.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = tas.diff, maxcell = Inf) +
  geom_sf(data = st_simplify(ecoregions, preserveTopology = TRUE, dTolerance = 1000), color = gray(.5), fill = NA, size = 0.01) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +  coord_sf(xlim = c(box.sam[1], box.sam[3]), ylim = c(box.sam[2], box.sam[4]), expand = TRUE) +
  annotate(geom = "text", x = -74, y = 2, label = "Colombia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -72, y = -7.5, label = "Brazil", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -75, y = -5, label = "Peru", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -76.6, y = -0.75, label = "Ecuador", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -70.75, y = 7.4, label = "Venezuela", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -73.2, y = -3.5, label = "NAP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -77.6, y = -4.5, label = "CCP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -73.6, y = 11.5, label = "SMP", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = -70.9, y = 10.1, label = "CMP", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"), na.value = NA,
                       limits = c(min(values(tas.diff, na.rm = TRUE)),
                                  round(max(values(tas.diff, na.rm = TRUE)), digits = 1))) +
  scale_x_continuous(breaks = seq(-80, -70, 2),
                     labels = c(-80, -78, -76, -74, -72, -70)) + 
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = c(-10, -5, 0, 5, 10)) + 
  labs(x = NULL, y = NULL, title = "Neotropic") +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = paste0("\u0394", "T", " (\u00b0C)"),
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

tas.diff.sam.fig

sam.x.range <- ggplot_build(pr.diff.sam.fig)$layout$panel_params[[1]]$x_range
sam.y.range <- ggplot_build(pr.diff.sam.fig)$layout$panel_params[[1]]$y_range

tas.diff.sam.fig <- tas.diff.sam.fig + annotation_custom(ggplotGrob(inset), xmin = sam.x.range[1], xmax = sam.x.range[1]+4, ymax = sam.y.range[2], ymin = sam.y.range[2]-2)

tas.diff.sam.fig

#Standard error
tas.diff.se.sam.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = tas.diff.se, maxcell = Inf) +
  geom_sf(data = st_simplify(ecoregions, preserveTopology = TRUE, dTolerance = 1000), color = gray(.5), fill = NA, size = 0.01) +
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
  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"), na.value = NA,
                       limits = c(0, 0.5)) +
  scale_x_continuous(breaks = seq(-80, -70, 2),
                     labels = c(-80, -78, -76, -74, -72, -70)) + 
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = c(-10, -5, 0, 5, 10)) + 
  labs(x = NULL, y = "Temperature"
       ) +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = "\u00b0C",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

tas.diff.se.sam.fig

# tas.diff.se.sam.fig <- tas.diff.se.sam.plot + annotation_custom(ggplotGrob(sam.inset), xmin = sam.x.range[1], xmax = sam.x.range[1]+4, ymax = sam.y.range[2], ymin = sam.y.range[2]-2)
# tas.diff.se.sam.fig


# Africa ------------------------------------------------------------------

#mean
tas.diff.afr.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = tas.diff, maxcell = Inf) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +  geom_sf(data = ecoregions, color = gray(.5), fill = NA, size = 0.01) +
  annotate(geom = "text", x = 39, y = 0, label = "Kenya", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35, y = 8, label = "Ethiopia", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35.5, y = -2.7, label = "Tanzania", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 36.5, y = 0.8, label = "EAMM", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = 39.5, y = 5.9, label = "EMM", size = 2.5, color = "black", fontface = "bold") +
  coord_sf(xlim = c(box.afr[1], box.afr[3]), ylim = c(box.afr[2], box.afr[4]), expand = TRUE) +
  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"), na.value = NA,
                       limits = c(min(values(tas.diff, na.rm = TRUE)),
                                  round(max(values(tas.diff, na.rm = TRUE)), digits = 1))) +
  scale_x_continuous(breaks = seq(35, 40, 2),
                     labels = c(35, 37, 39)) + 
  scale_y_continuous(breaks = seq(0, 10, 5),
                     labels = c(0, 5, 10)) + 
  labs(x = NULL, y = NULL, title = "Afrotropic") +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = paste0("\u0394", "T", " (\u00b0C)"),
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

tas.diff.afr.fig

# tas.diff.afr.fig <- tas.diff.afr.plot + annotation_custom(ggplotGrob(afr.inset), xmin = afr.x.range[1], xmax = afr.x.range[1]+3.31, ymax = afr.y.range[2], ymin = afr.y.range[2]-1.59)
# 
# tas.diff.afr.fig

#Standard error
tas.diff.se.afr.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = tas.diff.se, maxcell = Inf) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  geom_sf(data = ecoregions, color = gray(.5), fill = NA, size = 0.01) +
  annotate(geom = "text", x = 39, y = 0, label = "Kenya", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35, y = 8, label = "Ethiopia", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35.5, y = -2.7, label = "Tanzania", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 36.5, y = 0.8, label = "EAMM", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = 39.5, y = 5.9, label = "EMM", size = 2.5, color = "black", fontface = "bold") +
  coord_sf(xlim = c(box.afr[1], box.afr[3]), ylim = c(box.afr[2], box.afr[4]), expand = TRUE) +
  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"), na.value = NA,
                       limits = c(0, 0.5)) +
  scale_x_continuous(breaks = seq(35, 40, 2),
                     labels = c(35, 37, 39)) + 
  scale_y_continuous(breaks = seq(0, 10, 5),
                     labels = c(0, 5, 10)) + 
  labs(x = NULL, y = NULL 
       ) +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = "\u00b0C",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

tas.diff.se.afr.fig

# tas.diff.se.afr.fig <- tas.diff.se.afr.plot + annotation_custom(ggplotGrob(afr.inset), xmin = afr.x.range[1], xmax = afr.x.range[1]+3.31, ymax = afr.y.range[2], ymin = afr.y.range[2]-1.59)
# tas.diff.se.afr.fig

# Papua New Guinea --------------------------------------------------------

#mean
tas.diff.r <- t(flip(stack(tas.diff), 2))
tas.diff.r.df <- as.data.frame(tas.diff.r, xy = TRUE) %>% drop_na()
names(tas.diff.r.df) <- c("x", "y", "X1")

tas.diff.png.fig <- ggplot() +
  geom_sf(data = world.rot, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_raster(aes(x = x, y = y, fill = X1), data = tas.diff.r.df) +
  geom_sf(data = crsg_rect.rot, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  geom_sf(data = crsg_rect.rot, fill = NA, color = "black", size = 2, linetype = "dashed") +
  geom_sf(data = eco.png.rot, color = gray(.5), fill = NA, size = 0.01) +
  annotate(geom = "text", x = -6.5, y = -140, label = "Indonesia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -7, y = -141.5, label = "Papua New Guinea", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -6.85, y = -148, label = "300 km", color = "black", size = 1.9) +
  annotate(geom = "text", x = -6, y = -135.8, label = "CRSG", size = 2.5, color = "black", fontface = "bold") +
  coord_sf(xlim = c(box.png[1], box.png[3]), ylim = c(box.png[2], box.png[4]), expand = TRUE) +
  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"),
                       limits = c(min(tas.diff.r.df$X1),
                                  round(max(tas.diff.r.df$X1), digits = 1))) +
  scale_x_continuous(breaks = seq(-8, -4, 2),
                     labels = c(-8, -6, -4)) + 
  scale_y_continuous(breaks = seq(-148, -136, 4),
                     labels = c(148, 144, 140, 136)) +
  labs(x = NULL, y = NULL, title = "Australasia") +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = paste0("\u0394", "T", " (\u00b0C)"),
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

tas.diff.png.fig

# tas.diff.png.fig <- tas.diff.png.plot + annotation_custom(ggplotGrob(png.inset), xmin = png.x.range[1], xmax = png.x.range[1]+2.3, ymax = png.y.range[2], ymin = -136.2)
# 
# tas.diff.png.fig

#standard error

tas.diff.se.r <- t(flip(stack(tas.diff.se), 2))
tas.diff.se.r.df <- as.data.frame(tas.diff.se.r, xy = TRUE) %>% drop_na()
names(tas.diff.se.r.df) <- c("x", "y", "X1")

tas.diff.se.png.fig <- ggplot() +
  geom_sf(data = world.rot, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_raster(aes(x = x, y = y, fill = X1), data = tas.diff.se.r.df) +
  geom_sf(data = crsg_rect.rot, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  geom_sf(data = eco.png.rot, color = gray(.5), fill = NA, size = 0.01) +
  annotate(geom = "text", x = -6.5, y = -140, label = "Indonesia", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -7, y = -141.5, label = "Papua New Guinea", size = 1.5, color = "grey22") +
  annotate(geom = "text", x = -6.85, y = -148, label = "300 km", color = "black", size = 1.9) +
  annotate(geom = "text", x = -6, y = -135.8, label = "CRSG", size = 2.5, color = "black", fontface = "bold") +
  coord_sf(xlim = c(box.png[1], box.png[3]), ylim = c(box.png[2], box.png[4]), expand = TRUE) +
  scale_fill_gradientn(colours = c("white", "yellow", "orange", "red", "darkred"),
                       limits = c(0, 0.5)) +
  scale_x_continuous(breaks = seq(-8, -4, 2),
                     labels = c(-8, -6, -4)) + 
  scale_y_continuous(breaks = seq(-148, -136, 4),
                     labels = c(148, 144, 140, 136)) +
  labs(x = NULL, y = NULL
       ) +
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
  guides(fill = guide_colourbar(barwidth = 0.6, barheight = 7.5,
                                ticks = FALSE,
                                ticks.colour = "black",
                                draw.ulim = FALSE, draw.llim = FALSE,
                                frame.colour = "black",
                                title = "\u00b0C",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

tas.diff.se.png.fig

# tas.diff.se.png.fig <- tas.diff.se.png.plot + annotation_custom(ggplotGrob(png.inset), xmin = png.x.range[1], xmax = png.x.range[1]+2.3, ymax = png.y.range[2], ymin = -136.2)
# tas.diff.se.png.fig

# Produce paper figures with patchwork -------------------------------------

#Ensemble mean change of tas and pr (figure 2)
r1 <- (tas.diff.sam.fig | tas.diff.afr.fig | tas.diff.png.fig) + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "right")

r2 <- (pr.diff.sam.fig | pr.diff.afr.fig | pr.diff.png.fig) + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "right")

plot <- r1 / r2 + plot_annotation(tag_levels = list(c("A", "", "", "B", "", ""))) & theme(plot.tag = element_text(face = "bold",                                                                                                         size = 12),
                                                                                          plot.tag.position = "topleft")
plot

ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE1_NEW.pdf", device = "pdf", dpi = 1000, width = 16, height = 18, units = "cm")
ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE1_NEW.svg", device = "svg", dpi = 1000, width = 16, height = 18, units = "cm")
ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE1_NEW.png", device = "png", dpi = 1000, width = 16, height = 18, units = "cm")


#Ensemble standard error of tas and pr
r1 <- (tas.diff.se.sam.fig + labs(x = NULL, y = paste("Standard Error of", paste0("\u0394", "T")))| tas.diff.se.afr.fig | tas.diff.se.png.fig) + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "right")

r2 <- (pr.diff.se.sam.fig + labs(x = NULL, y = paste("Standard Error of", paste0("\u0394", "P"))) | pr.diff.se.afr.fig | pr.diff.se.png.fig) + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "right")

plot <- r1 / r2 + plot_annotation(tag_levels = list(c("A", "", "", "B", "", ""))) & theme(plot.tag = element_text(face = "bold",
                                                                                                                  size = 12),
                                                                                          plot.tag.position = "topleft")
plot

ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE_S1.pdf", device = "pdf", dpi = 1000, width = 16, height = 18, units = "cm")
ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE_S1.svg", device = "svg", dpi = 1000, width = 16, height = 18, units = "cm")
ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE_S1.png", device = "png", dpi = 1000, width = 16, height = 18, units = "cm")

