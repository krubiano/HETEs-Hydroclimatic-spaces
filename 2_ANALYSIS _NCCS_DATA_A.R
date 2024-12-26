################################################################################
################################################################################
###The script takes daily rasters by model/variable/experiment and returns rasters
###aggregated by month for approximately five years, clipped for the global 
###terrestrial surface
################################################################################
################################################################################

library(raster)
library(sf)
library(stringr)
library(lubridate)
library(ncdf4)
library(rgdal)

rasterOptions(progress = "text")
rasterOptions(timer = TRUE)
rasterOptions(chunksize = 1e+10)
rasterOptions(maxmemory = Inf)
rasterOptions(memfrac = 0.9)

rm(list=ls())
gc()


# Prepare input data ------------------------------------------------------------

# Load and reclassify previously generated raster global terrestrial boundaries
boundary_r <- raster("GIS/WORLD/land.tif")
boundary_r <- reclassify(boundary_r,
                         rcl = matrix(nrow = 2,
                                      ncol = 2,
                                      byrow = TRUE,
                                      data = c(c(0, 1), c( NA, NA))))

# Create variables to iterate

models <- dir("0_DATA/NEX-GDDP-CMIP6")
mymods <- c("IPSL-CM6A-LR") #Modify vector with the models you're interested in
models <- models[models %in% mymods]

myvars <- c("tas") #Modify vector with the variables you're interested in

myexps <- c("ssp585") #Modify vector with the experiments you're interested in to filter them

# df_files <- as.data.frame(matrix(nrow = 86+65))
# names(df_files) <- "year"
# df_files$year <- c(seq(1950, 2014, 1), seq(2015, 2100, 1))
# df_files$index <- c(1:65, 1:86)
# myyears <- c(2056:2060) #Modify with the starting and ending years of the files to filter
# myfiles <- df_files$index[df_files$year %in% myyears]

output_folder <- "1_PREP_DATA"
ext <- extent(-180, 180,-60, 90)
check_point <- as.character(c(seq(1955, 2014, 5), 2014, seq(2020, 2100, 5))) # Control some check points to save in disk

if (!output_folder %in% dir()) {
  dir.create(output_folder)
}


# Iteratively processes daily rasters by chunks -----------------------------

for (mod in 1:length(models)) {
  
  print(paste0("-(", Sys.time(), ")", " Model: ", models[mod]))
  
  variables <- dir(file.path("DATA/NEX-GDDP-CMIP6", models[mod]))
  
  if (exists("myvars")) {
    variables <- variables[variables %in% myvars]
  }
  
  if (!models[mod] %in% dir(output_folder)) {
    dir.create(file.path(output_folder, models[mod]))
  }
  
  
  for (var in 1:length(var)) {
    
    print(paste0("--(", Sys.time(), ")", " Variable: ", variables[var]))
    
    experiments <- dir(file.path("DATA/NEX-GDDP-CMIP6", models[mod], variables[var]))
    
    if (exists("myexps")) {
      experiments <- experiments[experiments %in% myexps]
    }
    
    if (!variables[var] %in% dir(file.path(output_folder, models[mod]))) {
      dir.create(file.path(output_folder, models[mod], variables[var]))
    }
    
    
    for (exp in 1:length(experiments)) {
      
      monthly_var <- stack()

      print(paste0("---(", Sys.time(), ")", " Experiment: ", experiments[exp]))
      
      files <- dir(file.path("DATA/NEX-GDDP-CMIP6", models[mod], variables[var], experiments[exp]))
      
      if (exists("myfiles")) {
        files <- files[myfiles]
      }
      
      if (!experiments[exp] %in% dir(file.path(output_folder, models[mod], variables[var]))) {
        dir.create(file.path(output_folder, models[mod], variables[var], experiments[exp]))
      }
      
      for (file in 1:length(files)) {
        
        print(paste0("----(", Sys.time(), ")", " File: ", files[file]))
        
        sdata <- stack(file.path("DATA/NEX-GDDP-CMIP6", models[mod], variables[var], experiments[exp], files[file]))
        print(paste0("(", Sys.time(), ")", " Stack done!"))
        
        nam <- names(sdata)
        
        # Rotates raster to standard geographical coordinates
        beginCluster()
        data <- clusterR(sdata, rotate)
        endCluster()
        
        print(paste0("(", Sys.time(), ")", " Rotate done!"))
        
        rm(sdata)
        gc()
        
        data <- setExtent(data, ext, keepres = TRUE)
        names(data) <- nam
        
        values(data)[is.na(values(boundary_r))] <- NA
        
        print(paste0("(", Sys.time(), ")", " Crop done!"))
        
        dates <- as.Date(str_replace_all(names(data), "X", ""), "%Y.%m.%d")
        months <- format(dates, "%m")
        months_int <- as.integer(months)
        year <- unique(format(dates, "%Y"))
        
        # aggregates rasters by months
        # Applies conversion factors and transformations according to the variable
        if (variables[var] == "pr") {
          data <- data*86400
          month_var <- stackApply(data, indices = months_int, fun = sum)
        } else if (variables[var] == "tas") {
          month_var <- stackApply(data, indices = months_int, fun = mean)
          month_var <- month_var-273.15
        } else if (variables[var] == "tasmax") {
          month_var <- stackApply(data, indices = months_int, fun = mean)
          month_var <- month_var-273.15
        } else if (variables[var] == "tasmin") {
          month_var <- stackApply(data, indices = months_int, fun = mean)
          month_var <- month_var-273.15
        } else if (variables[var] == "huss") {
          month_var <- stackApply(data, indices = months_int, fun = mean)
        } else if (variables[var] == "hurs") {
          month_var <- stackApply(data, indices = months_int, fun = mean)
        } else if (variables[var] == "rsds") {
          month_var <- stackApply(data, indices = months_int, fun = mean)
        }
        
        names(month_var) <- paste(year, unique(months))
        
        # Stacked iteratively month aggregated rasters
        monthly_var <- stack(monthly_var, month_var)
        
        # If the loop reaches a check point the file is saved to disk
        if (year %in% check_point) {
          
          writeRaster(x = monthly_var,
                      filename = file.path(output_folder, models[mod], variables[var], experiments[exp], paste0("TEMP_", year, ".tif")),
                      datatype = "FLT4S",
                      format = "GTiff",
                      overwrite = TRUE)
          
          monthly_var <- stack()
          gc()
        }
        
        rm(data)
        rm(sdata)
        rm(month_var)
        gc()
        
        print(paste0("-----","(", Sys.time(), ")", " file: ", files[file], " has been processed!"))
        
      }
      
    }
    
  }
  
}
