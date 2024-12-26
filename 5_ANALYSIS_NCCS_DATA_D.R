################################################################################
################################################################################
###The script aggregates the annual pixel values by 30-year epochs for recent
###historical and near future for each model.
################################################################################
################################################################################

library(raster)
library(sf)
library(terra)
library(dplyr)

memory.limit(size=560000)
rasterOptions(progress = "text")
rasterOptions(timer = TRUE)
rasterOptions(chunksize = 1e+10)
rasterOptions(maxmemory = Inf)
rasterOptions(memfrac = 0.9)
rasterOptions(tmptime = 72)

rm(list=ls())
gc()


# Load and prepare input data ---------------------------------------------

input_dir <- "3_ANNUAL_STACKS"
output_dir <- "4_ANNUAL_STACKS_SUMMARY"
output_dir_df <- "4A_ANNUAL_STACKS_SUMMARY_DF"

grassland <- st_read("GIS/WWF_TERR_ECO/wwf_tropical_montane.shp")
grassland_dis <- grassland %>% group_by(ECO_NAME, ECO_ID, G200_REGIO) %>% summarize()


# Aggregate the raw variables into two 30-year periods means by averaging pixel values --------

models <- dir(input_dir)
variables <- c("hurs", "pr", "rsds", "tas", "tasmax", "tasmin")
experiments <- c("historical", "ssp585")

for (mod in 1:length(models)) {
  
  for (var in 1:length(variables)) {
    
    for (exp in 1:length(experiments)) {
      
      r <- rast(list.files(file.path(input_dir, models[mod], variables[var], experiments[exp]), pattern = "\\.RDS$|\\.grd$|\\.tif$", full.names = TRUE))
      r <- r[[(nlyr(r)-29):nlyr(r)]]
      
      if (!dir.exists(file.path(output_dir, models[mod], variables[var], experiments[exp]))) {
        dir.create(file.path(output_dir, models[mod], variables[var], experiments[exp]), recursive = TRUE)}
      
      r <- tapp(r, index = rep(1, 30), fun = "mean",
                filename = file.path(output_dir, models[mod], variables[var], experiments[exp], paste0(variables[var], "_mod_summary_", experiments[exp], ".tif")),
                filetype = "GTiff",
                datatype = "FLT4S",
                overwrite = TRUE,
                memfrac = 0.9)
      
      r <- stack(r)
      
      #Extract values and create df
      r.data <- raster::extract(r, grassland_dis,
                                weights = TRUE,
                                df = TRUE,
                                exact = TRUE,
                                small = TRUE,
                                normalizeWeights = TRUE,
                                cellnumbers = TRUE)
      
      if (!dir.exists(file.path(output_dir_df, models[mod], variables[var], experiments[exp]))) {
        dir.create(file.path(output_dir_df, models[mod], variables[var], experiments[exp]), recursive = TRUE)}
      
      write.csv(r.data, file.path(output_dir_df, models[mod], variables[var], experiments[exp], paste0(variables[var], "_mod_summary_", experiments[exp], ".csv")))
      
    }
    
  }
  
}



# Compute the difference among reference and future periods by pix --------

for (mod in 1:length(models)) {
  
  for (var in 1:length(variables)) {
    
    r.norm.h <- stack(list.files(file.path(output_dir, models[mod], variables[var], "historical"),
                                 full.names = TRUE, recursive = TRUE, pattern = "\\.tif$"))
    
    r.norm.f <- stack(list.files(file.path(output_dir, models[mod], variables[var], "ssp585"),
                                 full.names = TRUE, recursive = TRUE, pattern = "\\.tif$"))
    
    r.norm.dif <- r.norm.f - r.norm.h
    
    if (!dir.exists(file.path(output_dir, models[mod], variables[var], "diff"))) {
      dir.create(file.path(output_dir, models[mod], variables[var], "diff"), recursive = TRUE)}
    
    #Save diff rasters
    raster::writeRaster(r.norm.dif,
                        filename = file.path(output_dir, models[mod], variables[var], "diff", paste0(variables[var], "_mod_summary_", "diff", ".tif")),
                        format = "GTiff",
                        datatype = "FLT4S",
                        overwrite = TRUE)
    
    #Extract df from rasters
    data.dif <- raster::extract(r.norm.dif, grassland_dis,
                                weights = TRUE,
                                df = TRUE,
                                exact = TRUE,
                                small = TRUE,
                                normalizeWeights = TRUE,
                                cellnumbers = TRUE)
    
    if (!dir.exists(file.path(output_dir_df, models[mod], variables[var], "diff"))) {
      dir.create(file.path(output_dir_df, models[mod], variables[var], "diff"), recursive = TRUE)
    }
    
    write.csv(data.dif, file.path(output_dir_df, models[mod], variables[var], "diff", paste0(variables[var], "_mod_summary_", "diff", ".csv")))
    
    rm(r.norm.h, r.norm.f, r.norm.dif, data.dif)
    gc()
    
  }
  
}


# Summarize data for model/variable/ecoregion -----------------------------
#Create a df sumarizing the values

experiments <- c("historical", "ssp585", "diff")

std.error <- function(x) sd(x)/sqrt(length(x))

for (exp in 1:length(experiments)) {

  for (var in 1:length(variables)) {
    
    if (!dir.exists(file.path(output_dir_df, "ENSEMBLE_SUMMARY", variables[var], experiments[exp]))) {
      dir.create(file.path(output_dir_df, "ENSEMBLE_SUMMARY", variables[var], experiments[exp]), recursive = TRUE)
    }

    df.list <- list()

    for (mod in 1:length(models)) {

      df.mod <- read.csv(list.files(file.path(output_dir_df, models[mod], variables[var], experiments[exp]), full.names = TRUE))
      names(df.mod) <- c("X", "ID", "cell", "X1", "weight")
      df.mod <- df.mod %>%
        group_by(ID) %>%
        summarise(wmean = weighted.mean(x = X1, w = weight))
      names(df.mod) <- c("ID", models[mod])

      df.list[[mod]] <- df.mod
    }

    df.ensemble <- do.call(cbind, df.list)
    df.ensemble <- df.ensemble[, !duplicated(colnames(df.ensemble))]
    df.ensemble$mean <- rowMeans(df.ensemble[2:7])
    df.ensemble$standard_dev <- apply(df.ensemble[2:7], 1, sd)
    df.ensemble$standard_err <- apply(df.ensemble[2:7], 1, std.error)

    write.csv(df.ensemble, file.path(output_dir_df, "ENSEMBLE_SUMMARY", variables[var], experiments[exp], paste0(variables[var], "_ens_summary_", experiments[exp], ".csv")))

  }

}
