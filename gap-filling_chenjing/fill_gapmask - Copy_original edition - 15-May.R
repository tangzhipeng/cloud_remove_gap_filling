# this code is used to fill gaps in Landsat 7
# written by Toby, 2 March, 2018
require("XML")
require("raster")
require("rgdal")

# clean the working environment
rm(list=ls())

fill_points <- function(ind_arr,input){
  # Written by Toby
  # May 15 2018
  # ind_arr is the return data frame from the function fill_gapmask
  # input is another NDVI file of TM or ETM+
  
  target <- ind_arr[[3]]
  if (length(ind_arr[[1]]) !=0 | length(ind_arr[[3]]) !=0) {
    output_name <- paste(substr( target,1,nchar(target)-13 ), substr(input,14,22),"fil.tif", sep="")
    ra_tar <- raster(target)
    ra_input <- raster(input)
    
    matri_ta <- as.matrix(ra_tar)
    matri_in <- as.matrix(ra_input)
    # at first, we define the size of neighbrhood is rad=2
    rad <- 2
    # we define the estimated number of m_class=2 
    m_class <-1.2
    # we define the minimum number of similar pixels as m_similar=20
    m_similar <- 20
    RMSD <- array(NA,c(50,50))
    DIST <- array(NA,c(50,50))
    nb_fil <- 0 # it calculates the number of filled pixels
    nb_unfil <- 0
    x_unfil <- c() # it records the x- y- of unfilled pixels
    y_unfil <- c()
    #let's enlarge the extend of each file for 20 pixels
    en <- 20
    n_row <- nrow(matri_ta)
    n_col <- ncol(matri_ta)
    target <- enlarge_matrix(matri_ta,en)
    input <- enlarge_matrix(matri_in,en)
    
    sum_na <- 0
    sum_all <- 0
    
    rm(matri_ta)
    rm(matri_in)
 
    # let's judge if pixels are in cloud or not, and only keep the in-cloud ones
    for (k in 1:length(ind_arr[[1]])) {
      i <- ind_arr[[1]][k]
      j <- ind_arr[[2]][k]
      if (i>10) {
        # firstly, count NA numbers in input
        sum_na <- length(which(is.na(input[ (i-10):(i+10),(j-10):(j+10) ])))
        sum_all <- length(input[ (i-10):(i+10),(j-10):(j+10) ])     
      }
      if (sum_na <= sum_all/2){
        #Let's search for the similar pixels
        count_similar <- 0
        while( (count_similar < m_similar) & (rad <= 10) ) {
          # this is for calulating the number of similar pixels and RMSD
          ##consider if the growing neighbour size meets the boundary
          # situation 1
          temp1 <- input[ (i-rad):(i+rad),(j-rad):(j+rad) ]
          dim(temp1) <- c(length(temp1))
          threshold <- sd(temp1,na.rm=TRUE)*2/m_class #standard deviation
          
          for ( x in (i-rad):(i+rad) ) {
            for(y in (j-rad):(j+rad) ){
              RMSD[x-(i-rad)+1,y-(j-rad)+1] <- abs(input[x,y]-input[i,j]) #RMSD countains all the 
              DIST[x-(i-rad)+1,y-(j-rad)+1] <- sqrt((x-i)^2+(y-j)^2)
            }
          }
          RMSD_nocent <- RMSD[1:(2*rad+1),1:(2*rad+1)]
          RMSD_nocent <- RMSD_nocent[-2*rad^2+2*rad+1]
          similar_RMSD <- RMSD_nocent
          similar_RMSD[similar_RMSD > threshold] <- NA #keep similar pixels
          
          DIST_nocent <- DIST[1:(2*rad+1),1:(2*rad+1)]
          DIST_nocent <- DIST_nocent[-2*rad^2+2*rad+1]
          similar_DIST <- DIST_nocent + similar_RMSD - similar_RMSD
          
          common_inp <- input[ (i-rad):(i+rad),(j-rad):(j+rad) ]
          common_inp <- common_inp[-2*rad^2+2*rad+1]
          similar_inp <- common_inp[!is.na(similar_RMSD)]
          
          common_tar <- target[(i-rad):(i+rad),(j-rad):(j+rad)]
          common_tar <- common_tar[-2*rad^2+2*rad+1]
          similar_tar <- common_tar[!is.na(similar_RMSD)]
          
          #similar <- 
          count_similar <- length(which(!is.na(similar_RMSD)))
          rad <- rad +1
        }
        rad <- rad-1
        
        #normolnize weight
        similar_weig <- similar_RMSD * similar_DIST
        similar_weig <- similar_weig[!is.na(similar_RMSD)] 
        similar_weig[is.na(similar_tar)]=NA
        similar_weig <- similar_weig / sum(similar_weig, na.rm = TRUE)
        
        
        #keep similar pixels and remove NA value
        
        #first, calculate L1 in target
        if ( length(which(!is.na(similar_tar)))==0 ){# find no similar pixels
          #reset rad
          rad <- 6
          temp2 <- target[ (i-rad):(i+rad),(j-rad):(j+rad) ]
          dim(temp2) <- c(length(temp2))
          mean_tar <- mean(temp2,na.rm=TRUE)
          standev_tar <-sd(temp2,na.rm=TRUE)
          
          temp3 <- input[ (i-rad):(i+rad),(j-rad):(j+rad) ]
          dim(temp3) <- c(length(temp3))
          mean_inp <- mean(temp3,na.rm=TRUE)
          standev_inp <- sd(temp3,na.rm=TRUE)
          gain <- standev_tar/standev_inp
          bias <- mean_tar - mean_inp * gain
          target[i,j] <- gain * input[i,j]+bias
          
        }
        
        else {
          center1 <- sum( similar_weig * similar_tar, na.rm = TRUE )
          center2 <- input[i,j] + sum( similar_weig * (similar_tar-similar_inp), na.rm = TRUE )
          R1 <- sum( abs(similar_inp-input[i,j]), na.rm = TRUE )
          R2 <- sum( abs(similar_inp-similar_tar), na.rm = TRUE)
          T1 <- 1/R1/(1/R1+1/R2)
          T2 <- 1/R2/(1/R1+1/R2)
          target[i,j] <- T1*center1+T2*center2
        }
        
        #reset rad to the original value
        rad <- 2
        RMSD <- array(NA,c(50,50))
        DIST <- array(NA,c(50,50))
        if (!is.na(target[i,j])) {
          nb_fil <- nb_fil + 1
        }
        else {
          x_unfil[ nb_unfil+1 ] <- i
          y_unfil[ nb_unfil+1 ] <- j
          nb_unfil <- nb_unfil + 1
        }
        rm(RMSD_nocent)
        rm(DIST_nocent)
        
      }
      else {
        x_unfil[ nb_unfil+1 ] <- i
        y_unfil[ nb_unfil+1 ] <- j
        nb_unfil <- nb_unfil + 1
      }
    }
    target_gapfilled <- raster(target[(en+1):(n_row+en), (en+1):(n_col+en)])
    extent(target_gapfilled) <- extent(ra_tar)
    projection(target_gapfilled) <- projection(ra_tar)
    writeRaster(target_gapfilled, output_name ,"GTiff", datatype='FLT4S')
    print(paste0("this process filled the number of pixes: ",nb_fil))
    print(paste0("The remian unfilled pixes are: ",nb_unfil))
    result <- list(x_unfil, y_unfil, output_name)
    return(result)
  }
  else {
    print("The unfilled array is empty and we did nothing")
    result <- list(c(), c(), c())
    return(result)
  }
}

fill_gapmask <- function(gap_mask_file,target,input){ 
  # Written by Toby
  # May 15 2018
  # gap_mask is the name of gap mask where values greater than 0 is the mask we want.
  # Before running this code, we need to clip gap_mask, target, and input files to make them in the same size
  # target is the NDVI file which needs to fill.
  # input is another NDVI file of TM or ETM+
  
  enlarge_matrix <- function(matri,en){
    #this function is used to enlarge the matrix by a size of "en"
    n_row <- nrow(matri)
    n_col <- ncol(matri)
    blank_matri <- array( NA,c(n_row+2*en,n_col+2*en) )
    blank_matri[1:n_row, 1:n_col] <- matri #topleft
    blank_matri[1:n_row, (2*en+1):(n_col+2*en)] <- matri #topright
    blank_matri[(2*en+1):(n_row+2*en), 1:n_col] <- matri #bottomleft
    blank_matri[(2*en+1):(n_row+2*en), (2*en+1):(n_col+2*en)] <- matri #bottomright
    blank_matri[(en+1):(n_row+en), (en+1):(n_col+en)] <- matri
    return(blank_matri)
  }
  
  startTime <- Sys.time()
  cat("Start time", format(startTime),"\n")
 
  output_name <- paste(substr(target,1,16), substr(input,14,22),"1st.tif", sep="")
 
  ra_gapma <- raster(gap_mask_file)
  ra_tar <- raster(target)
  ra_input <- raster(input)
  
  matri_ga <- as.matrix(ra_gapma)
  matri_ta <- as.matrix(ra_tar)
  matri_in <- as.matrix(ra_input)
  #check if these three files are in the same row and col
  
  # at first, we define the size of neighbrhood is rad=2
  rad <- 2
  # we define the estimated number of m_class=2 
  m_class <-1.2
  
  # we define the minimum number of similar pixels as m_similar=20
  m_similar <- 20
  RMSD <- array(NA,c(50,50))
  DIST <- array(NA,c(50,50))
  nb_fil <- 0 # it calculates the number of filled pixels
  nb_unfil <- 0
  x_unfil <- c() # it records the x- y- of unfilled pixels
  y_unfil <- c()
  
  
  #let's enlarge the extend of each file for 20 pixels
  en <- 20
  n_row <- nrow(matri_ga)
  n_col <- ncol(matri_ga)
  
  matri_gapma <- enlarge_matrix(matri_ga,en)
  target <- enlarge_matrix(matri_ta,en)
  input<- enlarge_matrix(matri_in,en)


  rm(matri_ga)
  rm(matri_ta)
  rm(matri_in)
  
  
  # judge in gap_mask and target files
  # circle in gapmask
  for ( j in (en+1):(n_col+en) ) {
    i <- en+1
    while (i <= n_row+en) {
      
      if ( (matri_gapma[i,j] > 0) & 
                        !is.na(matri_gapma[i,j]) ){
        up <- i
        #search the upper and lower line of maskgap. upper boundary is matri_gapma[up,j]
        #unsolved: if i==0
        
        while(matri_gapma[i,j] > 0 & 
              !is.na(matri_gapma[i,j]) &
              i <= n_row+en) 
          { i <- i+1 } # lower boundary is matri_gapma[dn,j]
        dn <- i-1 # lower boundary is matri_gapma[dn,j]
        
        # let's handle the pixels between up and down
        # NB! we should consider the boundry of up and dn
        # let's judge if pixels are in cloud or not, and only keep the in-cloud ones
        if ( !is.element(NA,target[ (up-6):(up-1),j ]) | 
             !is.element(NA,target[ (dn+1):(dn+6),j ]) ) {
          
          # do with the in-cloud pixels
          # let's judge if the 10-pixel (maximun) neighbourhood
          # between up and dn are NA in input pixels
          
          # firstly, count NA numbers in input
          sum_na <- length(which(is.na(input[ (up-10):(dn+10),(j-10):(j+10) ])))
          sum_all <- length(input[ (up-10):(dn+10),(j-10):(j+10) ])     
          if (sum_na < sum_all/2){          
            for (t in up:dn) {
              #Let's search for the similar pixels
              count_similar <- 0

              while( (count_similar < m_similar) & (rad <= 10) ) {
                # this is for calulating the number of similar pixels and RMSD
                ##consider if the growing neighbour size meets the boundary
                # situation 1
                temp1 <- input[ (t-rad):(t+rad),(j-rad):(j+rad) ]
                dim(temp1) <- c(length(temp1))
                threshold <- sd(temp1,na.rm=TRUE)*2/m_class #standard deviation
                
                for ( x in (t-rad):(t+rad) ) {
                  for(y in (j-rad):(j+rad) ){
                    RMSD[x-(t-rad)+1,y-(j-rad)+1] <- abs(input[x,y]-input[t,j]) #RMSD countains all the 
                    DIST[x-(t-rad)+1,y-(j-rad)+1] <- sqrt((x-t)^2+(y-j)^2)
                  }
                }
                RMSD_nocent <- RMSD[1:(2*rad+1),1:(2*rad+1)]
                RMSD_nocent <- RMSD_nocent[-2*rad^2+2*rad+1]
                similar_RMSD <- RMSD_nocent
                similar_RMSD[similar_RMSD > threshold] <- NA #keep similar pixels
                
                DIST_nocent <- DIST[1:(2*rad+1),1:(2*rad+1)]
                DIST_nocent <- DIST_nocent[-2*rad^2+2*rad+1]
                similar_DIST <- DIST_nocent + similar_RMSD - similar_RMSD
                
                common_inp <- input[ (t-rad):(t+rad),(j-rad):(j+rad) ]
                common_inp <- common_inp[-2*rad^2+2*rad+1]
                similar_inp <- common_inp[!is.na(similar_RMSD)]
                
                common_tar <- target[(t-rad):(t+rad),(j-rad):(j+rad)]
                common_tar <- common_tar[-2*rad^2+2*rad+1]
                similar_tar <- common_tar[!is.na(similar_RMSD)]
                
                #similar <- 
                count_similar <- length(which(!is.na(similar_RMSD)))
                rad <- rad +1
              }
              
              rad <- rad-1
              
              #normolnize weight
              similar_weig <- similar_RMSD * similar_DIST
              similar_weig <- similar_weig[!is.na(similar_RMSD)] 
              similar_weig[is.na(similar_tar)]=NA
              similar_weig <- similar_weig / sum(similar_weig, na.rm = TRUE)
              
              #keep similar pixels and remove NA value
              
              #first, calculate L1 in target
              if ( length(which(!is.na(similar_tar)))==0 ){# find no similar pixels
                #reset rad
                rad <- 6
                temp2 <- target[ (t-rad):(t+rad),(j-rad):(j+rad) ]
                dim(temp2) <- c(length(temp2))
                mean_tar <- mean(temp2,na.rm=TRUE)
                standev_tar <-sd(temp2,na.rm=TRUE)
                
                temp3 <- input[ (t-rad):(t+rad),(j-rad):(j+rad) ]
                dim(temp3) <- c(length(temp3))
                mean_inp <- mean(temp3,na.rm=TRUE)
                standev_inp <- sd(temp3,na.rm=TRUE)
                gain <- standev_tar/standev_inp
                bias <- mean_tar - mean_inp * gain
                target[t,j] <- gain * input[t,j]+bias
            
              }
              
              else {
                center1 <- sum( similar_weig * similar_tar, na.rm = TRUE )
                center2 <- input[t,j] + sum( similar_weig * (similar_tar-similar_inp), na.rm = TRUE )
                R1 <- sum( abs(similar_inp-input[i,j]), na.rm = TRUE )
                R2 <- sum( abs(similar_inp-similar_tar), na.rm = TRUE)
                #R1 <- sum( abs(similar_tar-input[t,j]), na.rm = TRUE )
                #R2 <- sum( abs(similar_tar-similar_inp), na.rm = TRUE)
                T1 <- 1/R1/(1/R1+1/R2)
                T2 <- 1/R2/(1/R1+1/R2)
                target[t,j] <- T1*center1+T2*center2
              }
              
              #reset rad to the original value
              rad <- 2
              RMSD <- array(NA,c(50,50))
              DIST <- array(NA,c(50,50))
              rm(RMSD_nocent)
              rm(DIST_nocent)
              if (!is.na(target[t,j])) {
                nb_fil <- nb_fil + 1
              }
              else {
                x_unfil[ nb_unfil+1 ] <- t
                y_unfil[ nb_unfil+1 ] <- j
                nb_unfil <- nb_unfil + 1
              }
            }

          }##
          else {
            x_unfil[ (nb_unfil+1) : (nb_unfil+1+dn-up) ] <- up:dn
            y_unfil[ (nb_unfil+1) : (nb_unfil+1+dn-up) ] <- j
            nb_unfil <- nb_unfil + dn-up+1
          }
        }
        
      }
      i <- i+1
    }
  }
  target_gapfilled <- raster(target[(en+1):(n_row+en), (en+1):(n_col+en)])
  extent(target_gapfilled) <- extent(ra_tar)
  projection(target_gapfilled) <- projection(ra_tar)
  writeRaster(target_gapfilled, output_name, "GTiff", datatype='FLT4S')
  print(paste0("this process filled the number of pixes: ",nb_fil))
  print(paste0("The remian unfilled pixes are: ",nb_unfil))
  
  if ( nb_unfil != 0 ){
    result <- list(x_unfil, y_unfil, output_name, nb_fil+nb_unfil)
    }
  else {
    result <- list(c(),c(),c(),nb_fil+nb_unfil)
    }
  
  return(result)
  timeDiff <- Sys.time() - startTime
  cat("\nProcessing time", format(timeDiff), "\n")
  #rm(matri)  # this is the final step in each circle as the matri takes lots of membery
}


startTime <- Sys.time()
cat("Start time", format(startTime),"\n")
Projectdirectory<- "E:/R_scripts/gap-filling_chenjing"  # Main project folder
setwd(Projectdirectory)
#files <- list.files()

####################################################  
#test for band3
gap_mask_files <- 'LE71670622011267_mask.tif'

target = 'LE71670622011267_band3_c1.tif'
input1 = 'LT51670622011211_band3_c1.tif'
input2 = 'LT51670622011307_band3_c1.tif'

index_unfil <- fill_gapmask(gap_mask_files,target,input1)
index_unfil2 <- fill_points(index_unfil,input2)
timeDiff <- Sys.time() - startTime
cat("\nProcessing time", format(timeDiff), "\n")

####################################
#test for band4
target = 'LE71670622011267_band4_c1.tif'
input1 = 'LT51670622011211_band4_c1.tif'
input2 = 'LT51670622011307_band4_c1.tif'

index_unfil_ <- fill_gapmask(gap_mask_files,target,input1)
index_unfil2_ <- fill_points(index_unfil_,input2)
timeDiff <- Sys.time() - startTime
cat("\nProcessing time", format(timeDiff), "\n")

####################################
#test for band5
target = 'LE71670622011267_band5_c1.tif'
input1 = 'LT51670622011211_band5_c1.tif'
input2 = 'LT51670622011307_band5_c1.tif'
index_unfil__ <- fill_gapmask(gap_mask_files,target,input1)
index_unfil2__ <- fill_points(index_unfil__,input2)
timeDiff <- Sys.time() - startTime
cat("\nProcessing time", format(timeDiff), "\n")

