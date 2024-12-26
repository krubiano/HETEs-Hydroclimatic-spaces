################################################################################
################################################################################
###The script joins the chunks of month aggregated rasters and produce a single
###raster file by y model/variable/experiment. Output rasters cointain 780 bands
###i.e. 12 bands/year
################################################################################
################################################################################

library(raster)
library(sf)
library(rgdal)
library(terra)

rasterOptions(progress = "text")
rasterOptions(timer = TRUE)
rasterOptions(chunksize = 1e+10)
rasterOptions(maxmemory = Inf)
rasterOptions(memfrac = 0.9)

terraOptions(progress = 2)
terraOptions(timer = TRUE)
terraOptions(memfrac = 0.9)

rm(list=ls())
gc()

input_folder <- "1_PREP_DATA"
output_folder <- "2_MONTHLY_STACKS"

models <- dir(input_folder)
mymods <- c("IPSL-CM6A-LR") #Modify vector with the models you"re interested in
models <- models[models %in% mymods]

myvars <- c("tas") #Modify vector with the variables you"re interested in

myexps <- c("ssp585") #Modify vector with the experiments you"re interested in

for (mod in 1:length(models)) {
  
  print(paste0("-(", Sys.time(), ")", " Model: ", models[mod]))
  
  variables <- dir(file.path(input_folder, models[mod]))
  
  if (exists("myvars")) {
    variables <- variables[variables %in% myvars]
  }
  
  for (var in 1:length(variables)) {
    
    print(paste0("--(", Sys.time(), ")", " Variable: ", variables[var]))
    
    experiments <- dir(file.path(input_folder, models[mod], variables[var]))
    
    if (exists("myexps")) {
      experiments <- experiments[experiments %in% myexps]
    }
    
    for (exp in 1:length(experiments)) {
      
      print(paste0("---(", Sys.time(), ")", " Experiment: ", experiments[exp]))
      
      if (!experiments[exp] %in% dir(file.path(output_folder, models[mod], variables[var]))) {
        dir.create(file.path(output_folder, models[mod], variables[var], experiments[exp]), recursive = TRUE)
      }
      
      files <- list.files(file.path(input_folder, models[mod], variables[var], experiments[exp]), full.names = TRUE, pattern = "\\.RDS$|\\.grd$|\\.tif$")
      
      for (file in 1:length(files)) {
        
        print(paste0("-----","(", Sys.time(), ")", " file: ", basename(files[file])))
        
        total <- 0
        
        if (substring(files[file], nchar(files[file]) - 3) == ".RDS") {
          temp_file <- rast(readRDS(files[file]))
        } else {
          temp_file <- rast(files[file])
        }
        
        n_occur <- data.frame(table(names(temp_file)))
        
        #Checks for errors or repeated months
        print(paste("Number of layers:", nlyr(temp_file)))
        print(paste("Number of years:", nlyr(temp_file)/12))
        print(paste("First file:", names(temp_file)[1]))
        print(paste("Last file:", names(temp_file)[nlyr(temp_file)]))
        print(paste("Repeated names:", names(temp_file)[names(temp_file) %in% n_occur$Var1[n_occur$Freq > 1]]))
        
        if ("raster_file" %in% ls()) {
          raster_file <- rast(list(raster_file, temp_file))
        } else {
          raster_file <- temp_file
        }
        
        total <- total + nlyr(temp_file)
        
        print(paste0("-----","(", Sys.time(), ")", " file: ", basename(files[file]), " has been processed!"))
        
      }
      
      print(paste0("-----","(", Sys.time(), ")", " Total nlayers: ", total))
      
      print(paste0("-----","(", Sys.time(), ")", " Exporting file: ", paste0(variables[var], "_month_", models[mod], "_", experiments[exp])))
      
      terra::writeRaster(x = raster_file,
                    filename = file.path(output_folder, models[mod], variables[var], experiments[exp], paste0(variables[var], "_month_", models[mod], "_", experiments[exp], ".tif")),
                    datatype = "FLT4S",
                    filetype = "GTiff",
                    overwrite = TRUE)
      
      print(paste0("-----","(", Sys.time(), ")", " File: ", paste0(variables[var], "_month_", models[mod], "_", experiments[exp]), " has been processed!"))
      
      rm(raster_file)
      rm(temp_file)
      gc()
    }
  }
}
