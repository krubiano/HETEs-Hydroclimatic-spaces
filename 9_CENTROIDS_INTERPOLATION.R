library(sf)
library(raster)

rm(list = ls())
dev.off()
gc()

input <- "4_ANNUAL_STACKS_SUMMARY"
output_df <- "4A_ANNUAL_STACKS_SUMMARY_DF"

# # Create centroids of polygons intersected by grid -----------------------

ecoregions <- as_Spatial(st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp"))
# ecoregions <- raster::shapefile("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp", encoding = "UTF-8")
ecoregions@data[which(ecoregions@data$ECO_NAME == "Santa Marta páramo"),][4] <- "Northern Andean páramo"

r <- raster(file.path(input, "ENSEMBLE", "biovar1", "historical", "biovar1_ens_1985-2014.tif"))
grid <- raster::rasterToPolygons(r, n = 4, na.rm = TRUE, dissolve = FALSE)
# writeOGR(obj = grid, dsn = "GIS/GRID/GRID.shp", layer = "GRID", driver = "ESRI Shapefile", overwrite_layer = TRUE)

eco.int <- raster::intersect(ecoregions, grid)

eco.int@data[which(eco.int@data$ECO_NAME == "Northern Andean pÃ¡ramo"),][4] <- "Northern Andean páramo"
eco.int@data[which(eco.int@data$ECO_NAME == "Cordillera de Merida pÃ¡ramo"),][4] <- "Cordillera de Merida páramo"
eco.int@data[which(eco.int@data$ECO_NAME == "Cordillera Central pÃ¡ramo"),][4] <- "Cordillera Central páramo"

eco.int <- st_cast(st_as_sf(eco.int), "POLYGON")
eco.int <- eco.int[c(1,4,11,22,23,24)]
eco.int$intersect_id <- seq(1, nrow(eco.int))
# st_write(eco.int, "GIS/GRID/ECO_POL_INTERSECT.shp", delete_layer = TRUE)

centroids <- sf::st_centroid(eco.int)
# st_write(centroids, "GIS/GRID/ECO_POL_CENTROIDS.shp", delete_layer = TRUE, driver = "ESRI Shapefile")

# # Extract interpolated data from raster to centroids ----------------------

# centroids <- st_read("GIS/GRID/ECO_POL_CENTROIDS.shp")
names(centroids) <- c("OBJECTID", "ECO_NAME", "G200_REGIO", "polygon_id", "ras_val", "geometry", "intersect_id")

vars <- c("biovar1", "biovar4", "biovar5", "biovar12", "biovar14", "biovar15") #Modify according to variables of interest
models <- dir(input)

# var <- vars[1] #Test
for (var in vars) {
  
  # mod <- models[1]
  for (mod in models) {
    
    bvar.h <- raster::extract(raster(file.path("4_ANNUAL_STACKS_SUMMARY", mod, var, "historical", paste0(var, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "1985-2014.tif"))),
                              centroids,
                              method = "bilinear",
                              sp = TRUE,
                              na.rm = TRUE)
    
    if (!dir.exists(file.path(output_df, "INTERPOLATED", var, mod, "historical"))) {
      dir.create(file.path(output_df, "INTERPOLATED", var, mod, "historical"), recursive = TRUE)}
    
    write.csv(bvar.h@data, file.path(output_df, "INTERPOLATED", var, mod, "historical", paste0(var, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "cen_1985-2014.csv")), fileEncoding = "UTF-8")
    
    
    bvar.585 <- raster::extract(raster(file.path("4_ANNUAL_STACKS_SUMMARY", mod, var, "ssp585", paste0(var, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "2071-2100.tif"))),
                                centroids,
                                method = "bilinear",
                                sp = TRUE,
                                na.rm = TRUE)
    
    if (!dir.exists(file.path(output_df, "INTERPOLATED", var, mod, "ssp585"))) {
      dir.create(file.path(output_df, "INTERPOLATED", var, mod, "ssp585"), recursive = TRUE)}
    
    write.csv(bvar.585@data, file.path(output_df, "INTERPOLATED", var, mod, "ssp585", paste0(var, ifelse(mod == "ENSEMBLE", "_ens_", "_mod_"), "cen_2071-2100.csv")), fileEncoding = "UTF-8")
    
  }
}
