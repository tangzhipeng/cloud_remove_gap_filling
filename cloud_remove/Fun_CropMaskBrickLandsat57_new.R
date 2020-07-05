#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
#*#* LANDSAT 7 AND LANDSAT 5 IMAGE WITH XMLVALUE ####
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*

Fun_CropMaskBrickLandsat57_new <- function(xml) {
  
  startTime <- Sys.time()
  cat("Start time", format(startTime),"\n")
  
  parse <- xmlParse(xml)
  root = xmlRoot(parse)
  
  # Check them manually if the xml file is still the same...
  # To check them use following steps parse <- xmlParse(images)
  # root = xmlRoot(parse)  # This images will come from below, use i=1 in console
  # fn_b1;fn_b2;fn_b3;fn_b4;fn_b5;fn_b7;fn_cfmask
  
  fn_b1 <- dir(pattern=".*.sr_band1.tif*")     #blue
  fn_b2 <- dir(pattern=".*.sr_band2.tif*")     #green
  fn_b3 <- dir(pattern=".*.sr_band3.tif*")     #red
  fn_b4 <- dir(pattern=".*.sr_band4.tif*")     #nir
  fn_b5 <- dir(pattern=".*.sr_band5.tif*")     #swir1
  fn_b7 <- dir(pattern=".*.sr_band7.tif*")     #swir2
  fn_cfmask <- dir(pattern=".*.pixel_qa.tif*") #pixel_qa band (cloud mask)
  
  # https://landsat.usgs.gov/landsat-surface-reflectance-quality-assessment
  # If using _pixel_qa.tif (band), 
  # Fill	1
  # Clear	66, 130
  # Water	68, 132
  # Cloud Shadow 72, 136
  # Snow/Ice	80, 112, 144, 176
  # Cloud	96, 112, 160, 176, 224
  # Low cloud confidence	66, 68, 72, 80, 96, 112
  # Medium cloud confidence	130, 132, 136, 144, 160, 176
  # High cloud confidence	224
  
  # If using sr_cloud_qa.tif (band), 			
  # DDV	1, 9
  # Cloud	2, 34
  # Cloud Shadow4, 12, 20, 36, 52
  # Adjacent to cloud	8, 12, 24, 40, 56
  # Snow	16, 20, 24, 48, 52, 56
  # Water	32, 34, 36, 40, 48, 52, 56
  
  #cfmask[cfmask == 2| cfmask == 4 | cfmask == 8 | cfmask == 34 |cfmask == 255] <- NA 
  #cfmask[cfmask == 1| cfmask == 2 | cfmask == 3 | cfmask == 4 | cfmask == 255] <- NA #(old format)
  #0> dark dense vegetation, 5>land/water, 4>snow
  tot_len = length(fn_b1)+length(fn_b2)+length(fn_b3)+length(fn_b4)+
    length(fn_b5)+length(fn_b7)+length(fn_cfmask)
  
  if (tot_len < 7) {
    message("Error! The bands are missing.")
  }
  
  else {
    metadata <- list.files(, pattern="MTL.txt")
    mtl <- readIniFile(metadata)
    mtl <- mtl[,-1]
    mtl[,2] <- gsub('\"', "", mtl[,2])
    imagename <- as.character(substr(as.character(mtl[mtl[,1]=="LANDSAT_SCENE_ID",2]),1,16 ))
    
    b1 <- raster(fn_b1)
    b2 <- raster(fn_b2)
    b3 <- raster(fn_b3)
    b4 <- raster(fn_b4)
    b5 <- raster(fn_b5)
    b7 <- raster(fn_b7)
    cfmask <- raster(fn_cfmask)
    # b1;b2;b3;b4;b5;b6;b7
    rm(fn_b1,fn_b2,fn_b3,fn_b4,fn_b5,fn_b7)
    
    
    # If using _pixel_qa.tif (band),
    cfmask[cfmask != 66 & cfmask != 68 & cfmask != 80 ] <- NA
    
    # 66 Clear, low confidence cloud
    # 68 Water, low confidence cloud  
    # 80 Snow/ice, low-confidence cloud
    
    # If using sr_cloud_qa.tif (band), 
    ###cfmask[cfmask != 1 & cfmask != 16 & cfmask != 32 & cfmask != 48 ] <- NA 
    # 1 Dark dense vegetation
    # 16 Snow 
    # 32 Water 
    # 48 Snow, water
    
    b1 <- mask(b1, cfmask)
    b2 <- mask(b2, cfmask)
    b3 <- mask(b3, cfmask)
    b4 <- mask(b4, cfmask)
    b5 <- mask(b5, cfmask)
    b7 <- mask(b7, cfmask)
    
    # b1 <- mask(b1, mask = studyarea, inverse = FALSE) # plot(b1)
    
    # plotRGB(image, r=7, g=5, b=3, stretch="lin")
    date <- gsub("-", "_", xmlValue(root[[1]][[4]][[1]]))
    
    # Lets see if some pixels are missed still in the cloud masking they are still in our image
    # So change those pixels (normally 20000 for clouds) to 0
    # Bad lines ( fill value) have value -9999, lets change them to 0
    # min value (-2000), lets change them to 0
    
    b1[b1==-9999| b1== 20000| b1== -2000]<- NA
    b2[b2==-9999| b2== 20000| b2== -2000]<- NA
    b3[b3==-9999| b3== 20000| b3== -2000]<- NA
    b4[b4==-9999| b4== 20000| b4== -2000]<- NA
    b5[b5==-9999| b5== 20000| b5== -2000]<- NA
    b7[b7==-9999| b7== 20000| b7== -2000]<- NA
    
    b1[b1< 0]<- NA # Valid range for reflectance is 0-10000 (scale factor 0.0001)
    b2[b2< 0]<- NA
    b3[b3< 0]<- NA
    b4[b4< 0]<- NA
    b5[b5< 0]<- NA
    b7[b7< 0]<- NA
    
    image2 <- brick(b1,b2,b3,b4,b5,b7)  ###
    fn_image2 <- paste(imagename,"_csmasked", ".tif", sep="")  ###
    #plotRGB(image2, r=6, g=5, b=3, stretch="lin")
    
    writeRaster(image2, fn_image2, "GTiff", datatype='FLT4S', options="INTERLEAVE=BAND", overwrite=TRUE)
    
    rm(image2, fn_image2)
    
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    #*#* Create NDVI ####
    #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*
    
    #NDVI<-  (b4-b3)/(b4+b3)  # Landsat 5 and 7 has band 4 as NIR and band 3 as red
    #writeRaster(NDVI, paste(imagename,"_ndvi_NoToCo",sep=""), "GTiff", overwrite=TRUE)
    
    timeDiff <- Sys.time() - startTime
    cat("\nProcessing time", format(timeDiff), "\n")
    
    rm(b1,b2,b3,b4,b5,b7)
  }
  
}
