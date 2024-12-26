library(sf)
library(dplyr)
library(terra)
library(raster)
library(fasterize)
library(rnaturalearth)
library(ggplot2)
library(tidyterra)
library(ggmap)
library(ggspatial)
library(RColorBrewer)
library(patchwork)

# Create intersection area raster -----------------------------------------

r <- raster(file.path("4_ANNUAL_STACKS_SUMMARY", "ENSEMBLE", "biovar1", "historical", "biovar1_ens_1985-2014.tif"))

ecoregions <- st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")
ecoregions[which(ecoregions$ECO_NAME == "Santa Marta páramo"),][4] <- "Northern Andean páramo"

grid <- st_read("GIS/GRID/GRID.shp")
grid$id <- seq(1:nrow(grid))

inters <- st_intersection(ecoregions, grid)

inters$int_area <- as.numeric(st_area(inters))

inters2 <- inters %>%
  group_by(id) %>% 
  summarize(int_area = sum(int_area))

grid_end <- left_join(grid, st_drop_geometry(inters2), by = "id")

grid_end$int_perc <- (grid_end$int_area/as.numeric(st_area(grid_end)) * 100)
  
raster_int <- fasterize(grid_end, r, field = "int_perc")

writeRaster(raster_int, filename = "GIS/GRID/ECO_GRID_INT.tif", format = "GTiff", overwrite = TRUE)


# Create plot -------------------------------------------------------------

# Load ecoregions shapefile and subset by location ------------------------

ecoregions <- st_read("/PASANTIA/GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")

r_int <- rast("GIS/GRID/ECO_GRID_INT.tif")

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


# South America -----------------------------------------------------------

int.sam.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = r_int, maxcell = Inf) +
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
  scale_fill_gradientn(colours = c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#006837", "#004529"), na.value = NA,
                       limits = c(0, 100)) +
  scale_x_continuous(breaks = seq(-80, -70, 2),
                     labels = c(-80, -78, -76, -74, -72, -70)) + 
  scale_y_continuous(breaks = seq(-10, 10, 5),
                     labels = c(-10, -5, 0, 5, 10)) + 
  labs(x = NULL, y = "Pixel area intersected by HETEs") +
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
                                title = "%",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        # panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.15),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

int.sam.fig


# Africa ------------------------------------------------------------------

int.afr.fig <- ggplot() +
  geom_sf(data = world, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_spatraster(data = r_int, maxcell = Inf) +
  geom_sf(data = ecoregions, color = gray(.5), fill = NA, size = 0.01) +
  geom_sf(data = rects, fill = NA, color = "black", linewidth = 0.35, linetype = "dashed") +
  coord_sf(xlim = c(box.afr[1], box.afr[3]), ylim = c(box.afr[2], box.afr[4]), expand = TRUE) +
  annotate(geom = "text", x = 39, y = 0, label = "Kenya", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35, y = 8, label = "Ethiopia", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 35.5, y = -2.7, label = "Tanzania", color = "grey22", size = 1.5) +
  annotate(geom = "text", x = 36.5, y = 0.8, label = "EAMM", size = 2.5, color = "black", fontface = "bold") +
  annotate(geom = "text", x = 39.5, y = 5.9, label = "EMM", size = 2.5, color = "black", fontface = "bold") +
  scale_fill_gradientn(colours = c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#006837", "#004529"), na.value = NA,
                       limits = c(0, 100)) +
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
                                title = "%",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        # panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.15),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

int.afr.fig



# Papua New Guinea --------------------------------------------------------

r_int.r <- t(flip(stack(r_int), 2))
r_int.r.df <- as.data.frame(r_int.r, xy = TRUE) %>% drop_na()
names(r_int.r.df) <- c("x", "y", "X1")

int.png.fig <- ggplot() +
  geom_sf(data = world.rot, color = gray(0.4), fill = gray(.93), size = 0.1) + 
  geom_raster(aes(x = x, y = y, fill = X1), data = r_int.r.df) +
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
  scale_fill_gradientn(colours = c("#FFFFE5", "#F7FCB9", "#D9F0A3", "#ADDD8E", "#78C679", "#41AB5D", "#238443", "#006837", "#004529"), na.value = NA,
                       limits = c(0, 100)) +
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
                         # pad_x = unit(2, "cm"),
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
                                title = "%",
                                title.vjust = 2)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        # panel.grid.major = element_line(color = gray(.5), linetype = "dashed", size = 0.15),
        axis.text = element_text(size = 6),
        legend.position = "none",
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 10, face = "bold"))

int.png.fig


# Produce paper figures with patchwork -------------------------------------

plot <- (int.sam.fig | int.afr.fig | int.png.fig) + plot_layout(nrow = 1, guides = "collect") & theme(legend.position = "right")
plot

ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE_INT.pdf", device = "pdf", dpi = 1000, width = 16, height = 9, units = "cm")
ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE_INT.svg", device = "svg", dpi = 1000, width = 16, height = 9, units = "cm")
ggsave(plot = plot, filename = "FIGURES/FINAL/FIGURE_INT.png", device = "png", dpi = 1000, width = 16, height = 9, units = "cm")
