rm(list=ls())
load("~/GitHub/SpatExtreme/HCDN_annual_max.RData")
library(evd)
library(ggplot2)
library(maps)
library(fields)
library(viridis)


Y       <- Y[,23:72]
year    <- year[23:72]
nmiss   <- rowSums(is.na(Y))
Y       <- Y[nmiss==0,]
s       <- s[nmiss==0,]
Y[Y<1]  <- 1
Y       <- log(Y)
X       <- (year-mean(year))/10
ns      <- nrow(Y)
nt      <- ncol(Y)

inits <- function(y,X,show=FALSE){
    library(ismev)
    fit    <- gev.fit(xdat=as.vector(y),ydat=X,mul=1,siglink=exp,show=show)
    return(fit)}

X <- matrix(X,50,1)

MLE <- matrix(0,ns,4)
for(i in 1:ns){
    fit     <- inits(Y[i,],X)
    MLE[i,] <- fit$mle
}

map <- function(s,y,main=""){
    dat <- data.frame(long=s[,1],lat=s[,2],MLE=y)
    g   <- ggplot(dat, aes(long, lat)) +
        borders("state") +
        geom_point(aes(colour = MLE)) +
        scale_colour_gradientn(colours = viridis(10)) +
        coord_map(projection = "albers", lat0 = 39, lat1 = 45) +
        xlab("")+ylab("")+labs(title=main)
    print(g)
}

map(s,MLE[,1],"GEV location - Intercept")
map(s,MLE[,2],"GEV location - Slope")
map(s,MLE[,3],"GEV log scale")
map(s,MLE[,4],"GEV shape")

dim(MLE)

Y1 = data.frame(MLE,lon = s$LONG_GAGE,lat = s$LAT_GAGE)
# Y1 = Y1[-drop,]
names(Y1)[1:4] = c('mu0','mu1','sigma','xi')
coordinates(Y1) = ~lon+lat
proj4string(Y1) = "+proj=longlat +datum=WGS84"
vario_mu0 = variogram(mu0~1,data = Y1)
vario_mu1 = variogram(mu1~1,data = Y1)
vario_sig = variogram(sigma~1,data = Y1)
vario_xi = variogram(xi~1,data = Y1)
plot(vario_mu0,main='mu')
plot(vario_mu1,main='mu')
plot(vario_sig,main='sigma')
plot(vario_xi,main='xi')


vario_mle = data.frame(distance = vario_mu0$dist, mu0 = vario_mu0$gamma, mu1 = vario_mu1$gamma, sig = vario_sig$gamma, xi = vario_xi$gamma)
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

