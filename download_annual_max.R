rm(list=ls())
#install.packages("dataRetrieval")
library(dataRetrieval)
library(stringr)

# The data come from
#   https://water.usgs.gov/osw/hcdn-2009/HCDN-2009_Station_Info.xlsx

# Download the data

 stn <- read.csv("HCDN/HCDN-2009_Station_Info.csv")

# Remove sites outside the CONUS

 outside <- stn[,6]%in%c(19,20,21)
 plot(stn[,8],stn[,7],col=1+outside)
 stn <- stn[!outside,]
 plot(stn[,8],stn[,7])

# Remove three sites with invalid(?) IDs

 badcode <- stn[,1]>100000000
 plot(stn[,8],stn[,7],col=1+badcode)
 stn <- stn[!badcode,]
 plot(stn[,8],stn[,7])

# Extract the important variables and discard the rest

 ID    <- as.character(stn[,1])
 ID    <- str_pad(ID, 8, pad = "0")
 s     <- stn[,8:7]
 HUC02 <- stn[,6]
 drain <- stn[,5]
 rm(stn,outside,badcode)

# Specify which data to download

 startDate   <- "1950-01-01"
 endDate     <- "2021-12-31"
 parameterCd <- c("00060")  # Discharge
 statCd      <- c("00003")  # Maximum
 nsites      <- length(ID)
 year        <- 1950:2020
 nyears      <- length(year)
 Y           <- matrix(NA,nsites,nyears)     
 min_obs_yr  <- 300

# Download the data

 for(i in 1:nsites){
   print(i)
   dat  <- readNWISdv(ID[i], parameterCd, startDate, endDate, statCd)
   dat  <- dat[!is.na(dat[,4]),]
   t    <- as.numeric(format(dat[,3], format = "%Y"))
   y    <- dat[,4]
   if(is.na(mean(y))){print("There are missing values")}
   for(j in 1:nyears){
     these <- which(t==year[j]) 
     if(length(these)>=min_obs_yr){
        Y[i,j] <- max(y[these])
     }
   } 
 } 

# Save the output
 
 rm(i,dat,y,j,t,these)
 rm(startDate,endDate,parameterCd,statCd,min_obs_yr)
 save.image("HCDN/HCDN_annual_max.RData")