library(beepr)
library(curl)

df <- read.csv("DOCUMENTOS/DOWNLOAD_FILES.csv")

#filter df if necessary
df <- subset(df, df$model == "MPI-ESM1-2-HR")
df <- subset(df, df$experiment == "ssp585")
df <- subset(df, df$variable == "tasmax")

output_root <- ("0_DATA/NEX-GDDP-CMIP6")

if (!unique(df$model) %in% dir(output_root)) {
  dir.create(file.path(output_root, unique(df$model)))
}

for (var in 1:length(unique(df$variable))) {
  
  if (!unique(df$variable)[var] %in% dir(file.path(output_root, unique(df$model)))) {
    dir.create(file.path(output_root, unique(df$model), unique(df$variable)[var]))
  }
  
  for (exp in 1:length(unique(df$experiment))) {
    
    if (!unique(df$experiment)[exp] %in% dir(file.path(output_root, unique(df$model), unique(df$variable)[var]))) {
      dir.create(file.path(output_root, unique(df$model), unique(df$variable)[var], unique(df$experiment)[exp]))
    }
    
    url_vector <- df$fileURL[df$variable == unique(df$variable)[var] & df$experiment == unique(df$experiment)[exp]]
    url_vector <- url_vector[72:86]
    
    for (link in 1:length(url_vector)) {
      
      filename <- basename(url_vector[link])
      
      if (filename %in% list.files(file.path(output_root,
                                             unique(df$model),
                                             unique(df$variable)[var],
                                             unique(df$experiment[exp])))) {
        next
      }
      
      destfile <- file.path(output_root,
                            unique(df$model),
                            unique(df$variable)[var],
                            unique(df$experiment)[exp],
                            filename)
      
      curl_download(url_vector[link], destfile, mode = "wb", quiet = FALSE)
      
      s_size <- as.numeric(httr::headers(httr::HEAD(url_vector[link]))[["Content-Length"]])
      
      d_size <- file.info(destfile)$size
      
      i <- 1
      while (d_size != s_size) {
        
        print(paste0("---", filename))
        print(paste("--File size:", s_size))
        print(paste("--Downloaded file size", d_size))
        
        beep(2)
        
        file.remove(destfile)
        
        curl_download(url_vector[link], destfile, mode = "wb", quiet = FALSE)
        
        i <- i+1
        
        if (i >= 3) {
          
          if (d_size != s_size) {
            
            beep(9)
            print(paste("-The file:", filename, "is corrupted"))
            
          }
          
          break
          
        }
        
      }
      
    }
    
  }

}
