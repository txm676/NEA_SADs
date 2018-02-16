
#load same packages as in grd_construction file


#####################################################################################
###########EXAMPLE SPATIAL REGRESSION CODE#######################################
###################################################################################

#run the grid and gambin functions
tes <- grid_sq_main(x, dis = 0.4, sens = FALSE, core = 8, SS = T, ss_V = 500, nr = 25) #run - 0.2 seems to be a good distance 

#remove NAs (grid squares with < 500 individuals)
tes <- tes[sapply(tes, function(x) all(!is.na(x)))] 

#for each returns a list; so turn into a df
matz <-  t(matrix(nrow = 9, unlist(tes))) 
matz <- as.data.frame(matz)
matz[, 1:3] <- apply(matz[, 1:3], 2, round, 2) #round
colnames(matz) <- c("Alpha", "N", "SR", "LAT", "LON", "Temp", "Seas", "Prec", "Chi")

##cut out outliers and bad x2 values
matz2 <- filter(matz, Alpha < 15)
matz2 <- filter(matz2, Chi >= 0.05)

#scale predictor variables
matz2 <- mutate(matz2, Temp = as.vector(scale(Temp)), Prec = as.vector(scale(Prec)),SR = as.vector(scale(SR)), N = as.vector(scale(N))) 


#get spatial data
cords <- cbind("x"= matz2$LON, "y"=matz2$LAT) #read in your x and y coord
kn4 <- spdep::knn2nb(spdep::knearneigh(cords, k=4, longlat = TRUE))
kn4_w <- spdep::nb2listw(kn4, style="W")

#run standard linear model and look at autocorrelation is residuals
mod <- lm(log(Alpha) ~ Temp +  Prec + SR, data = matz2)
spdep::moran.mc(residuals(mod), kn4_w, 999)#check autocorrelation of standard model


####################################################
#spatial regression model
###############################################################

#global model 
f1 <- 'log(Alpha) ~ Temp +  Prec + SR'
m1s <- spdep::lagsarlm(f1, data=matz2, kn4_w, tol.solve = 1e-30) #spatial lag model
m1e <- spdep::errorsarlm(f1, data=matz2, kn4_w, tol.solve=1.0e-30) #spatial error model

#compare spatial lag and spatial error models
AIC(m1s)
AIC(m1e)

#pick best and look at R2 and autocorrelation, and residual checks
spdep::moran.mc(residuals(m1e), kn4_w, 999)#re-check autocorrelation

#to get pseudo R2 and coefficient values
summary(m1e, Nagelkerke=TRUE) 

#look at residuals 
plot(m1e$fitted.values, residuals(m1e))#fitted vs residuals
resi <- (residuals(m1e))/  sqrt(m1e$s2)#residual var
qqnorm(resi)
qqline(resi)
hist(resi)
