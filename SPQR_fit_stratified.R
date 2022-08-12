rm(list = ls())
library(fields)
library(torch)
library(SPQR)
library(SpatialExtremes)
library(doParallel)
source("gen_data.R",local = T)
source("utils_MSP.R",local = T)
load("~/GitHub/SpatExtreme/HCDN_annual_max.RData")

Y1 = Y[,23:72]
dim(Y1)
howmanymissing = function(x){sum(is.na(x))}
missingnum = (apply(Y1,1,howmanymissing))
Y2 = Y1[missingnum==0,]
dim(Y2)
s1 = s#[HUC02==18,]
s2 = s1[missingnum==0,]
mins1 = min(s2[,1])
mins2 = min(s2[,2])
s2[,1] = s2[,1] - mins1
s2[,2] = s2[,2] - mins2
s_scale = max(s2)
s2 = s2/s_scale
dim(s2)
HUC02 = HUC02[missingnum==0]
s = s2#[HUC02==18,]
dim(s)
# test = data.frame(long = s[,1],lat = s[,2],huc = HUC02[HUC02==18])
# library(ggplot2)
# ggplot(test,aes(x=long,y=lat,col=huc,label=huc)) + geom_point() + geom_label()
# 
# maxdist = NA
# for(i in 1:489){
#     maxdist[i] = max(rdist(loc$nn_dist[,,i]))
#     print(paste(i,maxdist[i]))
# }
# if(cuda_is_available()){
#     use.GPU = T
# } else {
#     use.GPU = F
# }

n           <- nrow(s)
k           <- 15
# batch       <- 1
# fold        <- 7   # I run this seven times change this value from 1 to 7
# nfolds      <- 7
# which.fold <- rep(1:nfolds,n)[1:n]
# 
# print(batch)
# print(fold)

loc       <- create_locations(NS = 15, m = n, locs = s)

n_train   <- 200000
lr        <- 0.01
epochs    <- 50
n.knots   <- 15
n.hidden  <- c(30,20)

features <- function(delta,rho,r,Y){qnorm(rbind(delta,2*rho,r,Y))}


cl          <- makeCluster(10,outfile = 'Log.txt')
registerDoParallel(cl)
# registerDoSEQ()

fits <- foreach(id=2:489, .packages = c('SPQR','torch','SpatialExtremes')) %dopar% {
    mle.control     <- list(batch.size = 1000, epochs = epochs, lr = lr,
                            use.GPU = TRUE, early.stopping.epochs=20,save.name=paste0("mspgp_fit.pt"))
    tic     <- proc.time()
    min_delta <- 0.01
    max_delta <- 0.55
    set.seed(919*id + 1)
    nn        <- c(id,na.omit(loc$nn_id[id,]))
    dat       <- gen_data_evp_local3(n = n_train,s = as.matrix(loc$locs[nn,]),
                                     min_delta=min_delta,max_delta=max_delta)
    rho       <- dat$rho
    delta     <- dat$delta
    r = dat$r
    Y         <- dat$Y
    modelname <- paste('ffnn',id,'1',sep = '_')
    mle.control$save.name <- paste0("mspgp_fit_st",id,"batch 1.pt")
    y         <- Y[1,]
    X         <- t(features(delta=delta,rho=rho,r=r,Y=Y[-1,]))   
    fit       <- SPQR(X=X, Y=y, n.knots=n.knots, n.hidden=n.hidden,activation = "sigmoid",
                      method="MLE", control=mle.control, normalize=FALSE, verbose=F)
    save.SPQR(fit,modelname,path = 'EVP_HUC02')
    toc     <- proc.time()
    timetaken = (toc-tic)[3]
    print(paste('location',id,'batch 1 fitted in',timetaken,'seconds'))
    
    tic     <- proc.time()
    min_delta <- 0.45
    max_delta <- 0.99
    set.seed(919*id + 2)
    nn        <- c(id,na.omit(loc$nn_id[id,]))
    dat       <- gen_data_evp_local3(n = n_train,s = as.matrix(loc$locs[nn,]),
                                     min_delta=min_delta,max_delta=max_delta)
    rho       <- dat$rho
    delta     <- dat$delta
    r = dat$r
    Y         <- dat$Y
    modelname <- paste('ffnn',id,'2',sep = '_')
    mle.control$save.name <- paste0("mspgp_fit_st",id,"batch 2.pt")
    y         <- Y[1,]
    X         <- t(features(delta=delta,rho=rho,r=r,Y=Y[-1,]))
    fit       <- SPQR(X=X, Y=y, n.knots=n.knots, n.hidden=n.hidden,activation = "sigmoid",
                      method="MLE", control=mle.control, normalize=FALSE, verbose=F)
    save.SPQR(fit,modelname,path = 'EVP_HUC02')
    toc     <- proc.time()
    timetaken = (toc-tic)[3]
    print(paste('location',id,'batch 2 fitted in',timetaken,'seconds'))
}
stopCluster(cl)
