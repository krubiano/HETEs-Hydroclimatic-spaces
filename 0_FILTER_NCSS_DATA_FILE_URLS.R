library(dplyr)
library(stringr)

df <- read.csv("https://nex-gddp-cmip6.s3-us-west-2.amazonaws.com/gddp-cmip6-files.csv")

df$model <- ""
df$experiment <- ""
df$variable <- ""

fdf <- df %>%
  mutate(model = case_when(str_detect(fileURL, "ACCESS-ESM1-5") ~ "ACCESS-ESM1-5",
                           str_detect(fileURL, "EC-Earth3-Veg-LR") ~ "EC-Earth3-Veg-LR",
                           str_detect(fileURL, "EC-Earth3") ~ "EC-Earth3",
                           str_detect(fileURL, "HadGEM3-GC31-MM") ~ "HadGEM3-GC31-MM",
                           str_detect(fileURL, "CESM2") ~ "CESM2",
                           str_detect(fileURL, "IPSL-CM6A-LR") ~ "IPSL-CM6A-LR",
                           str_detect(fileURL, "MPI-ESM1-2-HR") ~ "MPI-ESM1-2-HR",
                           TRUE ~ ""))

ffdf <- fdf %>%
  mutate(experiment = case_when(str_detect(fileURL, "historical") ~ "historical",
                                str_detect(fileURL, "ssp126") ~ "ssp126",
                                str_detect(fileURL, "ssp585") ~ "ssp585",
                                TRUE ~ "")) 



fffdf <- ffdf %>%
  mutate(variable = case_when(str_detect(fileURL, "hurs") ~ "hurs",
                           str_detect(fileURL, "huss") ~ "huss",
                           str_detect(fileURL, "pr") ~ "pr",
                           str_detect(fileURL, "rlds") ~ "rlds",
                           str_detect(fileURL, "rsds") ~ "rsds",
                           str_detect(fileURL, "sfcWind") ~ "sfcWind",
                           str_detect(fileURL, "tas") ~ "tas",
                           str_detect(fileURL, "tasmax") ~ "tasmax",
                           str_detect(fileURL, "tasmin") ~ "tasmin",
                           TRUE ~ ""))

df_end <- fffdf[!(is.na(fffdf) | fffdf$model == "" | fffdf$experiment == "" | fffdf$variable == ""), ]
df_end <- na.omit(df_end)
df_end$fileURL <- trimws(df_end$fileURL, which = c("left"))

unique(df_end$model)

write.csv(df_end, "DOCUMENTOS/DOWNLOAD_FILES.csv")
