
library(ggplot2)
library(dplyr)
library(foreach)
library(doParallel)
library(cluster)
library(raster)
library(gambin)

##function 1: make the grid with a given distance between centre points
#x = a data frame of individual trees, which includes the plot id, the tree species,
#and the lat-lon co-ordinates of the plot (see the example "dat.csv" file for an example)
##sens = standard method (FALSE) or sensitivity analysis (TRUE) which moves top left corner
##dis = distance in latitudinal degrees
#as the length of a degree of longitude varies with latitude,
#this function accounts for this by varying the degrees of longtitude 
#used at each latitudinal band. 

#returns a list with three elements:
#1) gridP2 = a matrix of lat-lon co-ordiantes of the bottom left co-ordinates of the grid squares
#2)allSol = the longitude dis equivalent for each latitudinal band in the grid landscape
#3)latList = the latitude of bottom left corner of each grid square. Same as unique(gridP2$LAT)


grid_sq1 <- function(x, dis = 0.4, sens = FALSE){ ##sens = standard method (FALSE) or sensitivity analysis (TRUE) which moves top left corner
  
  ######################################
  #####make grid of fixed/random starting point######
  #########################################
  
  ######Find most extreme points from the empirical data to build the grid
  
  #South
  ws <- which.min(x$LAT)[1]
  xS <- filter(x, UN == x$UN[ws])#get unique id of this plot
  
  #North
  wn <- which.max(x$LAT)[1]
  xN <- filter(x, UN == x$UN[wn])#get unique id of this plot
  
  #east
  we <- which.max(x$LON)[1]
  xE <- filter(x, UN == x$UN[we])#get unique id of this plot
  
  #west
  ww <- which.min(x$LON)[1]
  xW <- filter(x, UN == x$UN[ww])#get unique id of this plot
  
  ##create grid matrix with lat-lon co-ordinates, with extreme points as boundaries
  ##create boundaries
  
  gridP <- matrix(ncol = 2, nrow = 0)
  colnames(gridP) <- c("LAT", "LON")
  
  #standard non-sensitivity analysis method
  if (sens == FALSE){
    
    #find top left corner (where N and W boundaries meet)
    TL <- c(xN$LAT[1], xW$LON[1])
    #then other corners
    TR <- c(xN$LAT[1], xE$LON[1])
    BR <- c(xS$LAT[1], xE$LON[1])
    BL <- c(xS$LAT[1], xW$LON[1])
  }
  
  #for sensitivity analysis that moves the grid starting points by random amount
  if (sens == TRUE){
    #random amount of variation to lat and lon of top left corner of up to 2.5
    ranV <- runif(1, 0.1, 2.5)
    cat(paste("random variation = ", ranV), "\n")
    
    #find top left corner (where N and W boundaries meet)
    TL <- c(xN$LAT[1], xW$LON[1])
    #now add and minus the random variation
    TL[1] <- TL[1] + ranV
    TL[2] <- TL[2] - ranV
    
    #then other corners (have checked all these for sense)
    TR <- c((xN$LAT[1] + ranV), xE$LON[1])
    BR <- c(xS$LAT[1], xE$LON[1])
    BL <- c(xS$LAT[1], (xW$LON[1] - ranV))
    
  }
  
  ######################################################
  #work out how much lon needs to change across latitude
  ######################################################
  
  #total height of the grid square landscape
  ble <- (ceiling(abs(TL[1] - BL[1]))) 
  #seg2 = the latitude of bottom left corner of each grid square
  seg2 <- seq(TL[1] - ble, (TL[1]), dis) 
  
  
  #function to convert degrees to radians
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  #work out distance of 1 degree of longtitude given the latitude. So one degree of lat "always" roughly 111km
  #but this function gives you the distance of one degree of longitude at a given lat band
  longDis <- function(lat) {111.320*cos(deg2rad(lat))}
  
  lat_dis <- dis * 111 #height of grid sq based on lat degrees (dis) we choose; one degree of latitude is ~ 111km
  
  #int function 
  #returns a sequence of longitude values relating to the longitude of each grid square along a latitudinal band
  int_fun5 <- function(lat_dis, lat_band){
    #work out the length of longitude needed to get a distance of 111km (to match lat)
    sol_dis <- lat_dis / longDis(lat_band) #gives the distance in longitude degrees needed to give same km distance as for latitude (i.e. 22.2 km for 0.2 degrees)
    Tle <- (ceiling(abs(TL[2] - TR[2])))  #total width of the grid landscape; the +2 makes sure the grid goes beyond the last point 
    seg <- seq((TL[2]), TL[2] + Tle, sol_dis) 
    return(seg)}
  
  #just get the distance multiplier for each lat band
  int_fun6 <- function(lat_dis, lat_band){
    #work out the length of longitude needed to get a distance of 111km (to match lat)
    #gives the distance in longitude degrees needed to give same km distance as for latitude
    #if lat_band = 0 (i.e. equator) then sol_dis should equal dis.
    sol_dis <- lat_dis / longDis(lat_band) 
    return(sol_dis )}
  
  #apply int_fun6 to each grid square (bottom corner), i.e. lat band
  #this then gives the longtitudinal version of 'dis' for each lat band
  allSol <- vapply(seg2, int_fun6, lat_dis = lat_dis, FUN.VALUE = numeric(1))
  
  seg2 <- seg2[which(seg2 <= TL[1])]#crop if needed
  
  #for each grid square (bottom corner) i.e. lat band, get the longitude co-ordinates
  #of the grid square corners. 
  lonList <- lapply(seg2, int_fun5, lat_dis = lat_dis) #lon values for each lat band
  latList <- seg2
  
  if (length(latList) != length(allSol)) stop ("lat and lon vector lengths do not match")
  
  tot_val <- (lapply(lonList, length) %>% unlist() %>% sum()) #total number of co-ords
  
  gridP2 <- matrix(ncol = 2, nrow = tot_val)
  
  if (length(latList) != length(lonList)) stop ("lat-lon issue: A")
  
  #for loops to go across each lat band starting from the bottom, and then for each lat band,
  #go across all the lon values in the stored list and save the lat-lon coordiantes in a matrix
  
  k <- 1 #tracker
  
  for (i in seq_along(latList)){
    latD <- latList[i]
    
    for (j in 1:length(lonList[[i]])){
      gridP2[k ,] <- c(latD, lonList[[i]][j])
      k <- k + 1
    }
  }
  
  gridP2 <- as.data.frame(gridP2)
  colnames(gridP2) <- c("LAT", "LON")
  
  resList <- list("gridP2" = gridP2, "allSol" = allSol, "latList" = latList)
  
  return(resList)}


########################################################
####fit Gambin to all of these grids####################
########################################################

##this function fits the gambin model to each of the samples created using the
#information from the grid_sq1 function. Takes the starting lat-lon co-ordinates
#and then builds the grid square using the 'dis' value provided (default = 0.4). dis refers
#the distance in latitude (i.e. the north-south boundary of a grid square). the distance in longitude
#(i.e. the east-west boundary of a grid square) varies across latitudinal bands and the value
#for each latitudinal band is taken from allSol.
#It runs a check to see if the sample has at least 500 individuals
#and if not returns NA.
#arguments are as for grid_sq1, except for 'core' which is the number of cores
#to use in the parallel processing, 'SS' which is whether or not to use subsampling
#standardisation when fitting the gambin model (ss_V is the number of individuals to subsample,
#is set to 500 in the main study. 'nr' is the number of times to repeat the subsampling before the mean
#alpha value is taken

#the function returns for each grid square either an NA (if less than 500 individuals in the square)
#or the alpha value, climatic variables, and number of individuals and species for the pooled data from
#that square

grid_sq2 <- function(gridP2 = gridP2, allSol = allSol, latList = latList,
                     x = x, dis = 0.4, core = 8, SS = F, ss_V = NULL, nr = 50){
    ##parallel processing code
    cores = core
    cl = makeCluster(cores); on.exit(stopCluster(cl))
    registerDoParallel(cl)
    i = 1 #Dummy line for RStudio warnings
    
    ##main parallel for loop
    matz = foreach(i=seq(from=1, to=nrow(gridP2), by=1))  %dopar% { 
      library(gambin)
      library(dplyr)
      
      #get the lat and lon of the ith plot
      lat1 <- gridP2$LAT[i]
      lon1 <- gridP2$LON[i]
      
      #work out the latitudinal band and then extract the correct longitudinal 'dis' for that lat band
      wh1 <- which(latList == lat1)
      l1 <- allSol[wh1]
      
      
      #make your grid square using the lat-lon co-ordinates of the bottom corner (lat1, lon1) and adding to
      #these the distance values for latitude (dis) and longitude (l1); and filter out all monitoring plots in this square
      dum <- filter(x, LAT >= lat1 & LAT < (lat1 + dis) & LON < (lon1 + l1) & LON >= lon1);
      #remember that in NEA longitude range is negative (~ -67 at extreme east to -96 at extreme west points)
      
      if (nrow(dum) < 500) return(NA) #sample size condition (has to return rather than 'next' as foreach)
      
      #get climate data (mean for all plots within local grid);
      #nb some sites close to coasts have no climate data in world clim
      #so return na: so filter these out
      temp1 <- ifelse(anyNA(dum$Temp), mean(na.omit(dum$Temp)), mean(dum$Temp))
      seas1 <- ifelse(anyNA(dum$Seas), mean(na.omit(dum$Seas)), mean(dum$Seas))
      pre1 <- ifelse(anyNA(dum$Prec), mean(na.omit(dum$Prec)), mean(dum$Prec))
      
      #create the abundance dataset & get alpha: subsample if necessary
      ab <- table(dum$SPCD) %>% unlist() %>% as.vector()
      if(SS == F){
        alp <- suppressMessages(gambin::fit_abundances(ab))
        alph <- round(alp$alpha, 2)
        chi <- summary(alp)$ChiSq$p.value
      }
      if(SS == T) {
        alp <- replicate(nr, summary(suppressMessages(gambin::fit_abundances(ab, subsample = ss_V)))) 
        aa <- c()
        chch <- c()
        for (i in 1:ncol(alp)){
          aa[i] <- alp[,i]$alp
          chch[i] <- alp[,i]$ChiSq$p.value
        }
        alph <- mean(aa) %>% round(2)
        chi <- mean(chch) %>% round(2)
      }
      
      #results vector
      c("alpha" = alph, "N" = nrow(dum), "SR" = length(ab), "LAT" = lat1, "LON" = lon1, "Temp" = temp1,
        "Seas" = seas1, "Prec" = pre1, "Chi" = chi) 
      
    }#eo i
  
    
    return(matz)
  }
  


##function combining 1 and 2

grid_sq_main <- function(x = x, dis = 0.4, sens = FALSE, core = 8, SS = F, ss_V = NULL, nr = 50){
  gg1 <- grid_sq1(x, dis, sens)
  gg2 <- grid_sq2(gridP2 = gg1[[1]], allSol = gg1[[2]], latList = gg1[[3]], x, dis, core, SS, ss_V, nr)
  return(gg2)
}

##example use
#test <- grid_sq_main(x,  dis = 0.4, sens = FALSE, core = 8, SS = T, ss_V = 500, nr = 10)

#############################################################
#####Multimodal Model Code##############################
##########################################################

#Similar to the above functions, but here both the unimodal and bimodal gambin models
#are fitted and compared using AIC. An option for subsampling is available but is not used
#in the main paper analyses as the alpha values are not compared here between samples. The 
#function arguments are the same as defined above. Also uses a sample size condition (NI) which
#by default is set to 500; any samples with fewer individuals than this return an NA. 

#Returns either an NA (for samples with < 500 individuals) or the AIC values and AIC weights for both
#the unimodal and bimodal gambin fits to a sample. The climate variables and number of species etc are also
#returned. 

grid_mm <- function(x = x, dis = 0.2, sens = FALSE, core = 4, NI = 500, SS = F, nr = 50){ #NI = no. indiv threshold
  gg <- grid_sq1(x, dis, sens)
  gg1 <- gg[[1]]
  allSol <- gg[[2]]
  latList <- gg[[3]]
  
  ##parallel processing code
  cores = core
  cl = makeCluster(cores); on.exit(stopCluster(cl))
  registerDoParallel(cl)
  i = 1 #Dummy line for RStudio warnings
  
  ##main parallel for loop
  matz = foreach(i=seq(from=1, to=nrow(gg1), by=1))  %dopar% { 
    library(gambin)
    library(dplyr)
    
    #get the lat and lon of the ith plot
    lat1 <- gg1$LAT[i]
    lon1 <- gg1$LON[i]
    
    #work out which lon dist to use
    wh1 <- which(latList == lat1)
    l1 <- allSol[wh1]
    
    
    #make your square using a set dis from the lat/lon point of the ith plot; and filter out all sites in square
    dum <- filter(x, LAT >= lat1 & LAT < (lat1 + dis) & LON < (lon1 + l1) & LON >= lon1) 
    
    if(nrow(dum) < NI) return(NA) #sample size condition (has to return rather than 'next' as foreach)
    
    #get climate data (mean for all plots within local grid);
    #nb some sites close to coasts have no climate data in world clim
    #so return na: so filter these out
    temp1 <- ifelse(anyNA(dum$Temp), mean(na.omit(dum$Temp)), mean(dum$Temp))
    seas1 <- ifelse(anyNA(dum$Seas), mean(na.omit(dum$Seas)), mean(dum$Seas))
    pre1 <- ifelse(anyNA(dum$Prec), mean(na.omit(dum$Prec)), mean(dum$Prec))
    
    #create the abundance dataset & get alpha
    
    if(SS == T){
      int3 <- function(dum, NI){
        abun <- table(dum$SPCD) %>% unlist() %>% as.vector()
        ab <- gambin:::create_octaves(abun, NI)
        
        fits1 = fit_abundances(ab, no_of_components = 1)#unimodal model
        fits2 = suppressMessages(fit_abundances(ab, no_of_components = 2, cores = 1))#bimodal model
        
        #calculate AIC, delta AIC and AIC weights for the two model fits
        aic <- c(AIC(fits1), AIC(fits2))
        
        aicb <- aic[which.max(aic)] - aic[which.min(aic)]
        
        if(which.min(aic) == 1) {daic <- c(0, aicb)} else{ daic <- c(aicb,0)}
        
        aicw1 <- vapply(daic, function(x) exp(-0.5 * x), numeric(1))
        
        aicw2 <- vapply(aicw1, function(x) x / (sum(aicw1)), numeric(1))
        
        #results vector
        res <- c("uniAIC" = aic[1], "biAIC" = aic[2], "uniW" = aicw2[1], "biW" = aicw2[2], "N" = nrow(dum), "SR" = length(abun), "LAT" = lat1, "LON" = lon1, "Temp" = temp1, "Prec" = pre1) 
        #remember N and SR are for the original sample, not the subsampled version (so will always be > 500 indivis etc)
        return(res)
      }#eo int3
      dfr <-  replicate(nr, int3(dum, NI))
      dfr2 <- dfr %>% as.data.frame() %>% t()
      
      #checks
      if(!any(apply(dfr2, 1, function(x) sum(x[3], x[4])) == 1)) warning("AIC weights don't sum to 1")
      if(anyNA(dfr)) warning("NAs in weights table")
      
      #return
      res <- colMeans(dfr2)
    }#eo if SS = T
    
    if(SS == F){
      ab <- table(dum$SPCD) %>% unlist() %>% as.vector() 
      fits1 = fit_abundances(ab, no_of_components = 1)
      fits2 = suppressMessages(fit_abundances(ab, no_of_components = 2, cores = 1))
      aic <- c(AIC(fits1), AIC(fits2))
      aicb <- aic[which.max(aic)] - aic[which.min(aic)]
      if(which.min(aic) == 1) {daic <- c(0, aicb)} else{ daic <- c(aicb,0)}
      aicw1 <- vapply(daic, function(x) exp(-0.5 * x), numeric(1))
      aicw2 <- vapply(aicw1, function(x)x / (sum(aicw1)), numeric(1))
      #results vector
      res <- c("uniAIC" = aic[1], "biAIC" = aic[2], "uniW" = aicw2[1], "biW" = aicw2[2], "N" = nrow(dum), "SR" = length(ab), "LAT" = lat1, "LON" = lon1, "Temp" = temp1, "Prec" = pre1) 
      #remember N and SR are for the original sample, not the subsampled version (so will always be > 500 indivis etc)
    }#eo ss=F
    
    res
    
  }#eo i
  return(matz)
}



