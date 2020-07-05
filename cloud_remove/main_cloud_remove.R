library("XML")
library("raster")
library("rgdal")
#library("maptools")
#library("gdata")
#library("doParallel")
# clear all the memory
rm(list=ls())

Maindirectory <- "E:/R_scripts/cloud_remove"
source(paste(Maindirectory, "Fun_CropMaskBrickLandsat57_new.R",sep="/"),echo=TRUE)
source(paste(Maindirectory, "Fun_CropMaskBrickLandsat8_new.R",sep="/"),echo=TRUE)

# open the dictionary where you store your Landsat files
subdic <- "E:/R_scripts/cloud_remove"
setwd(subdic)

# list all the dictionaries
files_LT5 <- list.files(,pattern = "LT05")
files_LE7 <- list.files(,pattern = "LE07")
files_LC8 <- list.files(,pattern = "LC08")


# process the Landsat
for (j in files_LC8) {
  
  # open the sub-dictionary
  setwd(file.path(subdic,j))
  
  # find the .xml file (Only one in each sub-dictionary)
  # remember do not open any .tif file in ArcGIS,otherwise there would be many .xml files
  file_xml <- list.files(,pattern = ".xml$")
  if (length(file_xml) != 0) {
    Fun_CropMaskBrickLandsat8_new(file_xml)
  }
  
}
