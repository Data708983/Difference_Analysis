rm(list = ls())
setwd('./')

if(!dir.exists("./Processed"))  dir.create("./Processed")
if(!dir.exists("./Processed/RawMatrix"))  dir.create("./Processed/RawMatrix")
if(!dir.exists("./Processed/RawData"))  dir.create("./Processed/RawData")
if(!dir.exists("./Processed/RawData/csv"))  dir.create("./Processed/RawData/csv")
if(!file.exists("./Data/gdc_download_20260315_132413.267760.tar.gz")){
  print("正在下载源数据\n")
  download.file("https://github.com/Data708983/Difference_Analysis/releases/download/source/gdc_download_20260315_132413.267760.tar.gz", 
                "./Data/gdc_download_20260315_132413.267760.tar.gz", mode = "wb")
}
tar_file = "./Data/gdc_download_20260315_132413.267760.tar.gz"
extract_dir = "./Processed/RawMatrix"
untar(tar_file, exdir = extract_dir)

#整理
library(data.table)
library(dplyr)
sample_sheet <- fread()