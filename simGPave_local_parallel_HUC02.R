rm(list=ls())
# library(geoR)
library(torch)
library(foreach)
library(doParallel)
library(SPQR)
library(SpatialExtremes)
load("~/GitHub/SpatExtreme/HCDN_annual_max.RData")
# cores = detectCores()
# cores

source("MCMC_veccia_approx_HUC02_SVC.R", local = T)
source("gradient.R", local = T)
source("gen_data.R", local = T)
source('utils_MSP.R', local = T)

Y       <- Y[,23:72]
s = s[,]
year    <- year[23:72]
nmiss   <- rowSums(is.na(Y))
Y       <- Y[nmiss==0,]
s       <- s[nmiss==0,]

s[,1]   <- s[,1]-min(s[,1])
s[,2]   <- s[,2]-min(s[,2])
s_scale <- max(s)
s       <- s/s_scale

Y[Y<1]  <- 1
Y       <- log(Y)
HUC02 = HUC02[nmiss==0]
# s = s[HUC02==18,]
# Y = Y[HUC02==18,]
X       <- (year-mean(year))/10
ns      <- nrow(Y)
nt      <- ncol(Y)


loc = create_locations(NS=15,m=ns,locs = s)
params_list_ensemble    <- list()
num_models              <- 2 # Number of models per location

features <- function(delta,rho,r,Y){qnorm(rbind(delta,2*rho,r,Y))}

for(q in 1:num_models){
    params_list <- vector(mode = 'list',length = ns)
    for(i in 2:ns){
        print(paste(q,i))
        objname         <- paste('ffnn',i,1,sep = '_')
        testmodel       <- load.SPQR(objname,path = 'EVP_HUC02_1model')
        ffnn_params     <- get.nn.params(testmodel)
        params_list[[i]]<- ffnn_params
    }
    params_list_ensemble[[q]] <- params_list
}

mean_deltas = rep(c(2,4,5,6,8)/10,2)
mean_deltas = c(0.2,0.8)
mean_deltas = qnorm(mean_deltas)
# rm(ffnn_params,params_list)
nsims=2
cl      <- makeCluster(2,outfile = 'Log.txt')
registerDoParallel(cl)
# registerDoSEQ()
output_app  <- foreach(sim = 1:nsims,.packages = c('Rfast','SpatialExtremes')) %dopar% {
    set.seed(sim)
    output_local  <- veccia_app_local(Y = Y,locs=loc,ffnn_params = params_list_ensemble,
                                      K=15,iters=20000,burn=1000,thin=1,update=100,mean_delta = mean_deltas[sim])
} 
stopCluster(cl)

par(mfrow=c(2,4))
for(i in 1:5){
    outHUC02=output_app[[i]]
    hist(outHUC02[[2]][,1],main = expression(delta),xlab = '')
    hist(outHUC02[[2]][,2],main = expression(rho),xlab = '')
    hist(outHUC02[[2]][,3],main = expression(r),xlab = '')
    hist(outHUC02[[2]][,4],main = expression(rho[2]),xlab = '')
    outHUC02=output_app[[i+1]]
    hist(outHUC02[[2]][,1],main = expression(delta),xlab = '')
    hist(outHUC02[[2]][,2],main = expression(rho),xlab = '')
    hist(outHUC02[[2]][,3],main = expression(r),xlab = '')
    hist(outHUC02[[2]][,4],main = expression(rho[2]),xlab = '')
}

for(i in 1:5){
    outHUC02=output_app[[i]]
    plot(1:20000,outHUC02[[2]][,1],main = expression(delta),xlab = '',type = 'l')
    plot(1:20000,outHUC02[[2]][,2],main = expression(rho),xlab = '',type = 'l')
    plot(1:20000,outHUC02[[2]][,3],main = expression(r),xlab = '',type = 'l')
    plot(1:20000,outHUC02[[2]][,4],main = expression(rho[2]),xlab = '',type = 'l')
    outHUC02=output_app[[i+1]]
    plot(1:20000,outHUC02[[2]][,1],main = expression(delta),xlab = '',type = 'l')
    plot(1:20000,outHUC02[[2]][,2],main = expression(rho),xlab = '',type = 'l')
    plot(1:20000,outHUC02[[2]][,3],main = expression(r),xlab = '',type = 'l')
    plot(1:20000,outHUC02[[2]][,4],main = expression(rho[2]),xlab = '',type = 'l')
}


par(mfrow=c(1,1))



mu0      <- 2*s[,1] #0.1444044
mu1     <- s[,2] #0.1858596
sig     <- 1
xi      <- 0.02
delta   <- .8
rho    <- .4
r = .6

# Generate data
truth <- c(mu0[21],mu1[21],sig,xi,delta,rho,r)

# length_common <-length(locs_common)
n       <- nrow(loc$locs)                               # total number of locations
s       <- loc$locs                                     # matrix of coordinates
k       <- ncol(loc$nn_id)                              # number of neighbors
K       <- length(ffnn_params[[length(ffnn_params)]])   # number of knots

nsims   <- 1
burn    <- 2000
iters   <- 5000
update  <- 500

Y = matrix(nrow = 34,ncol = 50)
# Generate datasets
for(i in 1:50){
    br_local <- t(rmaxstab(1, s, "brown", range = rhoW2rhoR(rho,1), smooth = 1))
    gev_mar <- matrix(rgev(nrow(s), 1, 1, 1), ncol = 1)
    tmp = cbind(r*br_local,(1-r)*gev_mar)
    R_s = apply(tmp,1,max)
    gp_local <- t(rgp(1, as.matrix(s), "powexp", nugget = r, range = rho, smooth = 1))
    br_exp   <- br_to_exp(R_s)
    gp_exp   <- gp_to_exp(gp_local,nugget = r)
    y        <- delta*br_exp + (1-delta)*gp_exp
    Y[,i]        <- phypoexp2(y, rate1 = 1/delta, rate2 = 1/(1-delta))
}
dim(Y)
years   <- 1972:2021
X_t     <- 0.1*(years-mean(years))
for(i in 1:34){
    location = mu0[i] + X_t*mu1[i]
    Y[i,] = qgev(Y[i,],location,sig,xi)
}
sigout=NA
xiout = NA
for(i in 1:34){
    sigout[i] = out$theta_gev[5000,i,3]
    xiout[i] = out$theta_gev[5000,i,4]
}
summary(sigout)
summary(xiout)
