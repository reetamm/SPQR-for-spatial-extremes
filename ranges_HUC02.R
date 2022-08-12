rm(list = ls())
library(gstat)
library(sp)
library(reshape2)
library(gdata)
library(fields)
library(dplyr)
library(maps)
library(ggplot2)
library(usmap)
load("~/GitHub/SpatExtreme/HCDN_annual_max.RData")


## US Map
MainStates <- map_data("state")
dim(Y)
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
# mins1 = min(s2[,1])
# mins2 = min(s2[,2])
# s2[,1] = s2[,1] - mins1
# s2[,2] = s2[,2] - mins2
# s_scale = max(s2)
# s2 = s2/s_scale
# s2
# Y=Y_log
s = s2
Y1 = apply(Y_log,1,quantile,0.9,na.rm=T)
length(Y1)
dim(s)
q95 = cbind(s,Y1)
names(q95) = c('lon','lat','streamflow')

ggplot() + geom_polygon(data=MainStates, aes(x=long, y=lat, group=group),color="black", fill="white" ) +
    geom_point(data = q95,aes(x=lon,y=lat,col=streamflow)) + coord_fixed(1.3) + #theme_bw() + 
    scale_color_gradient(low = "yellow", high = "darkblue")  +
    theme(legend.position = c(0.9,0.25),
          legend.title = element_blank(), legend.text = element_text(size = 14),
          axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank(),
          axis.title.x = element_blank(),axis.title.y = element_blank(),
          panel.border = element_blank(),panel.background = element_blank())


## chi 1 calculations
Y1 = Y_log
dim(Y1)
numpts = NA
for(i in 1:489)
{
    numpts[i] = sum(!is.na(Y1[i,]))
}
Y1_ranked = apply(Y1,1,rank,ties.method='average',na.last=T)
# View(Y1_ranked)
dim(Y1_ranked)
y1 = t(Y1_ranked)/(numpts+1)

u = c(0.75,0.85,0.95)
chi.array = array(dim=c(489,489,3))
for(i in 1:3){
    diag(chi.array[,,i]) = 1
}
for(i in 1:488)
    for(j in (i+1):489){
        avail = !is.na(Y1[i,]) & !is.na(Y1[j,])
        for(k in 1:3){
            if(max(y1[i,avail])>=1 | max(y1[j,avail])>=1)
                print(paste(i,j))
            numerator = sum(y1[i,avail]>u[k] & y1[j,avail]>u[k])
            denominator = sum(y1[j,avail]>u[k])
            if(denominator > 0){
                chi.array[i,j,k] = numerator/denominator
                chi.array[j,i,k] = numerator/denominator
            }
                
        }
    }
chi1 = chi.array[,,1]
chi2 = chi.array[,,2]
chi3 = chi.array[,,3]
chi1_v = as.vector(upperTriangle(chi1))
chi2_v = as.vector(upperTriangle(chi2))
chi3_v = as.vector(upperTriangle(chi3))
dist_s = as.vector(upperTriangle(rdist(s)))
dist_s = round(dist_s*111,0)
n=length(unique(dist_s))
chi1 = aggregate(chi1_v,list(dist_s),mean,na.rm=T)
chi2 = aggregate(chi2_v,list(dist_s),mean,na.rm=T)
chi3 = aggregate(chi3_v,list(dist_s),mean,na.rm=T)
names(chi1)  = c('dist','chi')
names(chi2)  = c('dist','chi')
names(chi3)  = c('dist','chi')

chi_df = rbind(chi1,chi2,chi3)
chi_df = as.data.frame(chi_df)
chi_df$u = rep(c('0.75','0.85','0.95'),each=n)
ggplot(chi_df,aes(x=dist,y=chi,col=u)) + geom_smooth(size=2) + theme_bw() + ylim(0,1)+
    labs(x='Distance (h)',y=expression(chi[h]),colour='Quantile (u)') +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )


#average variogram
ranges = NA
numyrs = NA
gamma = matrix(nrow = 50,ncol = 15)
Y1 = Y_log
Y1_ranked = apply(Y1,1,rank,ties.method='average',na.last=T)

for(i in 1:50){
    avail = which(!is.na(Y_log[,i]))
    numyrs[i] = length(avail)
    Y1 = Y1_ranked[i,avail]/50
    # Y1
    s1 = s[avail,]
    Y1_yr = data.frame(obs = Y1,lon = s1$LONG_GAGE,lat = s1$LAT_GAGE)
    coordinates(Y1_yr) = ~lon+lat
    proj4string(Y1_yr) = "+proj=longlat +datum=WGS84"
    vario_yr = variogram(obs~1,data = Y1_yr)
    gamma[i,] = vario_yr$gamma
    # print(plot(vario_yr))
    ranges[i] = fit.variogram(vario_yr, vgm(c("Exp","Sph","Mat")))$range[2]
}

dim(gamma)
gamma_mean = apply(gamma,2,mean)
plot(vario_yr$dist,gamma_mean,'b')
vario_df = data.frame(distance = vario_yr$dist, gamma = gamma_mean)
scaleFUN <- function(x) sprintf("%.2f", x)

ggplot(vario_df,aes(x=distance,y=gamma)) + geom_point(size=3) + geom_line() + theme_bw() +
    labs(x = 'Distance', y = 'Semivariance estimate') + ylim(0,0.1)+
    theme(legend.position = "bottom",
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    guides(color = guide_legend(override.aes = list(size = 2) ) ) 

vario_yr$gamma = gamma_mean
tmp = fit.variogram(vario_yr, vgm(c("Exp","Sph","Mat")))
plot(vario_yr,tmp)

plot(1950:2021,ranges)
summary(ranges[-14])
max(ranges[-14])/111 #6.7


## chi 2 plot
Y1 = Y_log
dim(Y1)
numpts = NA
s1 = rdist(s)*111
s1 = plyr::round_any(s1,10)
for(i in 1:489)
{
    numpts[i] = sum(!is.na(Y1[i,]))
}
Y1_ranked = apply(Y1,1,rank,ties.method='average',na.last=T)
# View(Y1_ranked)
dim(Y1_ranked)
y1 = t(Y1_ranked)/(numpts+1) 

u.seq = seq(0.5,0.99,by=.001)
n = length(u.seq)
dists =c(10,100,1000)
chi.seq = matrix(0,nrow = n,ncol = 3)
for(type in 1:3){
    count = 0
for(i in 1:488)
    for(j in (i+1):489){
        if(s1[i,j]==dists[type]){
            count=count+1
            avail = !is.na(Y1[i,]) & !is.na(Y1[j,])
            for(k in 1:n){
                u = u.seq[k]
                numerator = sum(y1[i,avail]>u & y1[j,avail]>u)
                denominator = sum(y1[j,avail]>u)
                if(denominator > 0){
                    chi.seq[k,type] =chi.seq[k,type] + numerator/denominator
                }
            }
        }
    }
    chi.seq[,type] = chi.seq[,type]/count
    }

head(chi.seq)
# par(mfrow=c(3,1))
# plot(u.seq,chi.seq[,1]) + title(main = 'dist = 0.1')
# plot(u.seq,chi.seq[,2]) + title(main = 'dist = 0.5')
# plot(u.seq,chi.seq[,3])+ title(main = 'dist = 1')
# par(mfrow=c(1,1))
chi = c(chi.seq[,1],chi.seq[,2],chi.seq[,3])
distsf = factor(dists)
levels(distsf) = c('10e0','10e1','10e2')
chi_df_2 = data.frame(chi=chi,u=rep(u.seq,3),dist=rep(distsf,each=n))
head(chi_df_2)
ggplot(chi_df_2,aes(x=u,y=chi,col=dist)) + geom_smooth(size=2) + theme_bw() +
    labs(x='Quantile (u)',y=expression(chi[u](h)),colour='Distance(h)') + ylim(0,1) +
    theme(legend.position = 'bottom',
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )

