load("~/GitHub/SPQR-for-spatial-extremes/MCMC_outputs.RData")
load("~/GitHub/SPQR-for-spatial-extremes/HCDN/HCDN_annual_max.RData")
library(fields)
library(torch)
library(SPQR)
library(SpatialExtremes)
library(doParallel)
source("gen_data.R",local = T)
source("utils_MSP.R",local = T)
post.sample <- seq(from=5001,to=20000,length.out=100)
post.sample <- round(post.sample,0)
dim(s)

trunc.data <- data_prep(Y,s)
Y       <- trunc.data$Y
s       <- trunc.data$s
nmiss   <- trunc.data$nmiss
ns      <- nrow(Y)
loc     <- create_locations(NS=15,m=ns,locs = s)

# Identify the sites in COlorado
stn     <- read.csv("HCDN/HCDN-2009_Station_Info.csv")
outside <- stn[,6]%in%c(19,20,21)
stn     <- stn[!outside,]
badcode <- stn[,1]>100000000
stn     <- stn[!badcode,]
state   <- stn[,9]
state   <- state[nmiss==0]
HUC02   <- HUC02[nmiss==0]
cluster <- which(HUC02==14 & state=='CO')
cluster <- unique(c(cluster,as.vector(loc$nn_id[cluster,])))

#Extract cluster parameters
dim(output_app[[1]][[2]])
theta_gev   <- output_app[[1]][[1]][post.sample,cluster,]
theta_sp    <- output_app[[1]][[2]][post.sample,]
theta_gev   <- abind::abind(theta_gev,output_app[[2]][[1]][post.sample,cluster,],along=1)
theta_sp    <- rbind(theta_sp,output_app[[2]][[2]][post.sample,])
y1          <- apply(Y[cluster,],1,quantile,0.9)

#Generate data for the cluster
exceedances <- matrix(NA,200,4)
for(i in 1:200){
    n=20000

    s_cl    <- as.matrix(s[cluster,])
    p       <- nrow(s_cl)
    rho     <- theta_sp[i,2]
    delta   <- theta_sp[i,1]
    r       <- theta_sp[i,3]
    R_s     <- matrix(nrow = n,ncol = p)
        
    br_local <- rmaxstab(n, s_cl, "brown", range = rhoW2rhoR(rho,1), smooth = 1)
    gev_mar  <- matrix(rgev(n*p, 1, 1, 1), ncol = p)
    for(locs in 1:p){
        tmp = cbind(r*br_local[,locs],(1-r)*gev_mar[,locs])    
        R_s[,locs]  <- apply(tmp,1,max)
        }
    gp_local<- rgp(n, s_cl, "powexp", nugget = r, range = rho, smooth = 1)
    br_exp  <- br_to_exp(R_s)
    gp_exp  <- gp_to_exp(gp_local,nugget = r)
    y       <- delta*br_exp + (1-delta)*gp_exp
    Y_cdf   <- phypoexp2(y[,1:7], rate1 = 1/delta, rate2 = 1/(1-delta))
    p       <- 7    
    n2      <- n/2
    Y1972   <- matrix(NA,n2,p)
    Y2021   <- matrix(NA,n2,p)
    
    for(j in 1:p){
        tmp1        <- Y_cdf[1:n2,j]
        tmp2        <- Y_cdf[(n2+1):n,j]
        Y1972[,j]   <- qgev(tmp1,loc = theta_gev[i,j,1]-2.45*theta_gev[i,j,2],scale = theta_gev[i,j,3], shape = theta_gev[i,j,4] )
        Y2021[,j]   <- qgev(tmp2,loc = theta_gev[i,j,1]+2.45*theta_gev[i,j,2],scale = theta_gev[i,j,3], shape = theta_gev[i,j,4] )
    }
    
    u   <- 0.5
    y1  <- y2 <- apply(t(Y[cluster[1:7],]),2,quantile,u)
    
    exceedances[i,1] <- sum(Y1972[,1]>y1[1] & Y1972[,2]>y1[2] & Y1972[,3]>y1[3] & Y1972[,4]>y1[4] & Y1972[,5]>y1[5] & Y1972[,6]>y1[6] & Y1972[,7]>y1[7])
    exceedances[i,2] <- sum(Y2021[,1]>y2[1] & Y2021[,2]>y2[2] & Y2021[,3]>y2[3] & Y2021[,4]>y2[4] & Y2021[,5]>y2[5] & Y2021[,6]>y2[6] & Y2021[,7]>y2[7])
    
    u   <- 0.9
    y1  <- y2 <- apply(t(Y[cluster[1:7],]),2,quantile,u)
    
    exceedances[i,3] <- sum(Y1972[,1]>y1[1] & Y1972[,2]>y1[2] & Y1972[,3]>y1[3] & Y1972[,4]>y1[4] & Y1972[,5]>y1[5] & Y1972[,6]>y1[6] & Y1972[,7]>y1[7])
    exceedances[i,4] <- sum(Y2021[,1]>y2[1] & Y2021[,2]>y2[2] & Y2021[,3]>y2[3] & Y2021[,4]>y2[4] & Y2021[,5]>y2[5] & Y2021[,6]>y2[6] & Y2021[,7]>y2[7])
    print(exceedances[i,])
}

exceedances <- exceedances/n2
hist(exceedances[,2]-exceedances[,1])
hist(exceedances[,4]-exceedances[,3])
sum(exceedances[,2]-exceedances[,1]>0,na.rm = T)/200

sd(exceedances[,4]-exceedances[,3]>0,na.rm = T)

stn2[cluster[1:7],]

mean(exceedances[,4]);sd(exceedances[,4])
