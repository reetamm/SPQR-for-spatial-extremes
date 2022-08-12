rm(list=ls())
load("HCDN/HCDN_annual_max.RData")
library(evd)
library(ggplot2)
library(maps)
library(fields)
library(viridis)
library(ismev)
library(sp)
library(gstat)

Y       <- Y[,23:72]            #choose the last 50 years
year    <- year[23:72]          #choose the last 50 years
nmiss   <- rowSums(is.na(Y))    #how many missing years at each locations?
Y       <- Y[nmiss==0,]         #choose locations with complete data
s       <- s[nmiss==0,]         #coordinates of locations with complete data
Y[Y<1]  <- 1                    #set minimum value to 1
Y       <- log(Y)               #log transform
X       <- (year-mean(year))/10 #scale the year vector
ns      <- nrow(Y)
nt      <- ncol(Y)

inits <- function(y,X,show=FALSE){
    fit    <- ismev::gev.fit(xdat=as.vector(y),ydat=X,mul=1,siglink=exp,show=show)
    return(fit)}

X <- matrix(X,50,1)

# Initial MLE estimates
MLE <- matrix(0,ns,4)
for(i in 1:ns){
    fit     <- inits(Y[i,],X)
    MLE[i,] <- fit$mle
}

Y1 = data.frame(MLE,lon = s$LONG_GAGE,lat = s$LAT_GAGE)
names(Y1)[1:04] <- c('mu0','mu1','sigma','xi')
coordinates(Y1) <- ~lon+lat
proj4string(Y1) <- "+proj=longlat +datum=WGS84"

# Variograms for the GEV MLEs
vario_mu0   <- variogram(mu0~1,data = Y1)
vario_mu1   <- variogram(mu1~1,data = Y1)
vario_sig   <- variogram(sigma~1,data = Y1)
vario_xi    <- variogram(xi~1,data = Y1)


# Data frame of MLE variograms
vario_mle = data.frame(distance = vario_mu0$dist, mu0 = vario_mu0$gamma, mu1 = vario_mu1$gamma, sig = vario_sig$gamma, xi = vario_xi$gamma)

# Plot variograms for MLE estimates
ggplot(vario_mle,aes(x=distance,y=mu0)) + geom_point(size=3) + geom_line()+ theme_bw() +
    labs(x = 'Distance', y = 'Semivariance estimate') + ylim(0,4)+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )

ggplot(vario_mle,aes(x=distance,y=mu1)) + geom_point(size=3) + geom_line()+ theme_bw() +
    labs(x = 'Distance', y = 'Semivariance estimate') + ylim(0,.02)+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )


ggplot(vario_mle,aes(x=distance,y=sig)) + geom_point(size=3) + geom_line() + theme_bw() +
    labs(x = 'Distance', y = 'Semivariance estimate') + ylim(0,.4)+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )


ggplot(vario_mle,aes(x=distance,y=xi)) + geom_point(size=3) + geom_line() + theme_bw() +
    labs(x = 'Distance', y = 'Semivariance estimate') + ylim(0,0.05)+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )