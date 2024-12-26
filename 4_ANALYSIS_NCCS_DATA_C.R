################################################################################
################################################################################
###The script aggregates monthly raster stacks into annual raster stacks and produce a single
###raster file by y model/variable/experiment. Output rasters cointain une band/year.
###output raster stacks are masked to Tropical Andean Ecosystem area. It also produces
###a DF by extracting pixel values mean by polygons.
################################################################################
################################################################################

library(terra)
library(dplyr)
library(raster)
library(sf)

rasterOptions(progress = "text")
rasterOptions(timer = TRUE)
rasterOptions(chunksize = 1e+10)
rasterOptions(maxmemory = Inf)
rasterOptions(memfrac = 0.9)
rasterOptions(tmptime = 72)

terraOptions(progress = 2)
terraOptions(timer = TRUE)
terraOptions(memfrac = 0.9)

rm(list=ls())
gc()

input <- "2_MONTHLY_STACKS"
output <- "3_ANNUAL_STACKS"
output_df <- "3A_ANNUAL_STACKS_DF"

models <- list.files(input)
mymods <- c("IPSL-CM6A-LR") #Modify vector with the models you"re interested in
models <- models[models %in% mymods]

myvars <- c("tas") #Modify vector with the variables you"re interested in

myexps <- c("ssp585")

grassland <- terra::vect("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")
sf.grassland <- st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")

for (mod in 1:length(models)) {
  
  variables <- list.files(file.path(input, models[mod]))
  
  if (exists("myvars")) {
    variables <- variables[variables %in% myvars]
  }
  
  for (var in 1:length(variables)) {
    
    experiments <- list.files(file.path(input, models[mod], variables[var]))
    
    if (exists("myexps")) {
      experiments <- experiments[experiments %in% myexps]
    }
    for (exp in 1:length(experiments)) {
      
      r <- terra::rast(list.files(file.path(input, models[mod], variables[var],
                                            experiments[exp]), full.names = TRUE,
                                  pattern = "\\.grd$|\\.tif$"))
      
      r.area <- terra::mask(r, grassland, touches = TRUE)
      
      rm(r)
      gc()
      
      if (!dir.exists(file.path(output, models[mod], variables[var], experiments[exp]))) {
        dir.create(file.path(output, models[mod], variables[var], experiments[exp]), recursive = TRUE)
      }
      
      if (variables[var] == "pr") {
        
        r.annual <- terra::tapp(r.area, index = rep(1:nlyr(r.area), each = 12),
                                fun = sum, cores = 7,  na.rm = TRUE,
                                filename = file.path(output, models[mod], variables[var],
                                                     experiments[exp],
                                                     paste0(variables[var], "_year_",
                                                            models[mod], "_",
                                                            experiments[exp],
                                                            ".tif")),
                                overwrite = TRUE,
                                wopt = list(filetype = "GTIff", datatype = "FLT4S"))
        
      } else {
        
        r.annual <- terra::tapp(r.area, index = rep(1:nlyr(r.area), each = 12),
                                fun = mean, cores = 7,  na.rm = TRUE,
                                filename = file.path(output, models[mod], variables[var],
                                                     experiments[exp],
                                                     paste0(variables[var], "_year_",
                                                            models[mod], "_",
                                                            experiments[exp],
                                                            ".tif")),
                                overwrite = TRUE,
                                wopt = list(filetype = "GTIff", datatype = "FLT4S"))
        
      }
      
      rm(r.area)
      gc()
      
      data.annual <- terra::extract(x = r.annual, y = grassland, fun = mean,
                                    na.rm = TRUE, list = FALSE, exact = TRUE)
      
      data.annual <- data.annual[,-1]
      
      df <- data.frame(x = unlist(as.data.frame(data.annual)))
      names(df) <- variables[var]
      
      df$polygon_id <- rep(1:134)
      
      if (experiments[exp] == "historical") {
        
        df$year <- rep(1950:2014)
        
      } else {
        
        df$year <- rep(2015:2100)
        
      }
      
      df$experiment <- rep(experiments[exp])
      
      data.annual <- merge(x = df, y = st_drop_geometry(sf.grassland[, c(22, 1, 8, 4, 11, 5)]), by = "polygon_id")
      
      if (!dir.exists(file.path(output_df, models[mod], variables[var], experiments[exp]))) {
        dir.create(file.path(output_df, models[mod], variables[var], experiments[exp]), recursive = TRUE)
      }
      
      write.csv(data.annual, file.path(output_df, models[mod], variables[var], experiments[exp],
                                       paste0(variables[var], "_", models[mod],
                                              "_", experiments[exp],".csv")))
      
      print(paste(models[mod], "-", variables[var], "-",
                  experiments[exp], " ready!"))
      
    }
  }
}