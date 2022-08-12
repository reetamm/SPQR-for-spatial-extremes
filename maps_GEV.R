library(ggplot2)
library(maps)
library(viridis)
library(usmap)
library(tidyr)
my.plots <- vector(2, mode='list')
length(output_app)

params = c('delta','rho','r','rho[2]')

par(mfrow=c(1,1))
outHUC02_1=output_app[[1]]
outHUC02_2=output_app[[2]]
plot(1:20000,outHUC02_1[[2]][,1],main = expression(delta),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='red',ylim = c(0.2,0.8))
lines(1:20000,outHUC02_2[[2]][,1],main = expression(delta),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='blue')
plot(1:20000,outHUC02_1[[2]][,2],main = expression(rho),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='red')
lines(1:20000,outHUC02_2[[2]][,2],main = expression(rho),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='blue')
plot(1:20000,outHUC02_1[[2]][,3],main = expression(r),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='red')
lines(1:20000,outHUC02_2[[2]][,3],main = expression(r),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='blue')
plot(1:20000,outHUC02_1[[2]][,4],main = expression(rho[2]),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='red')
lines(1:20000,outHUC02_2[[2]][,4],main = expression(rho[2]),ylab = 'Estimate',xlab = 'Iteration',type = 'l',col='blue')



par(mfrow=c(1,1))
out = output_app[[1]]
theta_sp = out$theta_sp[5001:20000,]
theta_gev = out$theta_gev[5001:20000,,]
out = output_app[[2]]
theta_sp = rbind(theta_sp,out$theta_sp[5001:20000,])
theta_gev = abind::abind(theta_gev,out$theta_gev[5001:20000,,],along = 1)
# theta_sp[,4] = pnorm(log(theta_sp[,4]))
theta_sp[,c(2,4)] = theta_sp[,c(2,4)]*(s_scale*111)
dim(theta_gev)
theta_gev[,,3] = log(theta_gev[,,3])
means_sp = apply(theta_sp,2,mean)
sd_sp = apply(theta_sp,2,sd)
means_gev = matrix(nrow = 489,ncol = 4)
for(i in 1:489){
    means_gev[i,] = apply(theta_gev[,i,], 2, mean)
}

MainStates <- map_data("state")
load("~/GitHub/SpatExtreme/HCDN_annual_max.RData")
Y1 = Y[,23:72]
dim(Y1)
howmanymissing = function(x){sum(is.na(x))}
missingnum = (apply(Y1,1,howmanymissing))
Y2 = Y1[missingnum==0,]
dim(Y2)
Y2[Y2<1] = 1
Y_log = log(Y2)
s1 = s#[1:29,]
s2 = s1[missingnum==0,]
means_gev = data.frame(means_gev,s2)
names(means_gev) = c('mu0','mu1','sigma','xi','lon','lat')

mu_exp = expression(mu[0])
ggplot() + geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),color="black", fill="white" ) +
    geom_point(data = means_gev,aes(x=lon,y=lat,col=mu0),size=2) + coord_fixed(1.3) + #theme_bw() + 
    scale_color_gradient(low = "yellow", high = "darkblue",name=mu_exp)  +
    theme(legend.position = c(0.92,0.25), legend.title.align=0.1,
          legend.text = element_text(size = 14),
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.border = element_blank(),panel.background = element_blank())

mu_exp = expression(mu[1])
ggplot() + geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),color="black", fill="white" ) +
    geom_point(data = means_gev,aes(x=lon,y=lat,col=mu1),size=2) + coord_fixed(1.3) + #theme_bw() + 
    scale_color_gradient(low = "yellow", high = "darkblue",name=mu_exp)  +
    theme(legend.position = c(0.92,0.25), legend.title.align=0.1,
          legend.text = element_text(size = 14),
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.border = element_blank(),panel.background = element_blank())
sig_exp = expression(sigma)
ggplot() + geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),color="black", fill="white" ) +
    geom_point(data = means_gev,aes(x=lon,y=lat,col=sigma),size=2) + coord_fixed(1.3) + #theme_bw() + 
    scale_color_gradient(low = "yellow", high = "darkblue",name=sig_exp)  +
    theme(legend.position = c(0.92,0.25), legend.title.align=0.1,
          legend.text = element_text(size = 14),
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.border = element_blank(),panel.background = element_blank())
xi1 = expression(xi)
ggplot() + geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),color="black", fill="white" ) +
    geom_point(data = means_gev,aes(x=lon,y=lat,col=xi),size=2) + coord_fixed(1.3) + #theme_bw() + 
    scale_color_gradient(low = "yellow", high = "darkblue",name=xi1)  +
    theme(legend.position = c(0.92,0.25), legend.title.align=0.1,
          legend.text = element_text(size = 14),
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.border = element_blank(),panel.background = element_blank())


dim(theta_gev)
mu1 = theta_gev[,,2]
mu1_positive = apply(mu1,2,function(x)mean(x>0))
summary(mu1_positive)
means_gev$mu1_positive = mu1_positive

prob = expression(P(mu[1]>0))
ggplot() + geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),color="black", fill="white" ) +
    geom_point(data = means_gev,aes(x=lon,y=lat,col=mu1_positive),size=2) + coord_fixed(1.3) + #theme_bw() + 
    scale_color_gradient(low = "yellow", high = "darkblue",name=prob)  +
    theme(legend.position = c(0.92,0.25), legend.title.align=-.1,
          legend.text = element_text(size = 14),
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.border = element_blank(),panel.background = element_blank())

dim(Y_log)
dim(means_gev)
gev_cdf = matrix(NA,489,50)
for(i in 1:50){
    loc = means_gev[,1] + X[i]*means_gev[,2]
    scale = (means_gev[,3])
    shape = means_gev[,4]
    gev_cdf[,i] = eva::pgev(Y_log[,i],loc,scale,shape)
}
for(i in 1:489)
    gev_cdf[i,] = sort(gev_cdf[i,])
quantiles = seq(0.01:0.99,by=0.02)
quantiles = quantile((1:50)/51,seq(0,1,length.out=50))
mean_cdf = apply(gev_cdf,2,mean)

gev_wide = data.frame(t(gev_cdf))
dim(gev_wide)
gev_wide$mean_cdf = mean_cdf
gev_wide$quantiles = quantiles
gev_long = gather(gev_wide,location,cdf,X1:X489)
names(gev_long)

p = ggplot(gev_long,aes(x=quantiles,y=cdf,fill=factor(location))) + geom_line(col='grey') + theme(legend.position = 'none')
p + geom_line(aes(y=mean_cdf,fill=factor(location))) + geom_abline(slope = 1,intercept = 0,color='red',linetype='dashed') +
    xlab('CDF') + ylab('Probability Integral Transform (PIT)')
summary(mu1_positive)
colors = rep('grey',489)
colors[mu1_positive<.21] = 'red'
colors[mu1_positive>.7] = 'blue'
par(mfrow=c(1,1))
plot(quantiles,gev_cdf[1,],col = colors[1],type = 'l')
for(i in 2:489)
    lines(quantiles,gev_cdf[i,],col = colors[i],type = 'l')

