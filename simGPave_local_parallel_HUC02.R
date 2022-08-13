rm(list=ls())
# library(geoR)
library(torch)
library(foreach)
library(doParallel)
library(SPQR)
library(SpatialExtremes)
load("HCDN/HCDN_annual_max.RData")
load("~/GitHub/SPQR-for-spatial-extremes/fitted_SPQR_models.RData")
# cores = detectCores()
# cores

source("MCMC_veccia_approx_HUC02_SVC.R", local = T)
source("gen_data.R", local = T)
source('utils_MSP.R', local = T)

trunc.data = data_prep(Y,s)
Y = trunc.data$Y
s = trunc.data$s

X       <- (year-mean(year))/10 # Standardize the effect of the year
ns      <- nrow(Y)
nt      <- ncol(Y)


loc <- create_locations(NS=15,m=ns,locs = s)
# params_list_ensemble    <- list()
# num_models              <- 2 # Number of models per location
# 
# features <- function(delta,rho,r,Y){qnorm(rbind(delta,2*rho,r,Y))}
# 
# for(q in 1:num_models){
#     params_list <- vector(mode = 'list',length = ns)
#     for(i in 2:ns){
#         print(paste(q,i))
#         objname         <- paste('ffnn',i,1,sep = '_')
#         testmodel       <- load.SPQR(objname,path = 'C:/Users/rmajumd3/Documents/GitHub/SpatExtreme/EVP_HUC02_1model')
#         ffnn_params     <- get.nn.params(testmodel)
#         params_list[[i]]<- ffnn_params
#     }
#     params_list_ensemble[[q]] <- params_list
# }

mean_deltas <- c(0.2,0.8) # 2 chains with 2 different starting values of delta
mean_deltas <- qnorm(mean_deltas)
nsims       <- length(mean_deltas)
cl          <- makeCluster(2,outfile = 'Log.txt')
registerDoParallel(cl)
# registerDoSEQ()
output_app  <- foreach(sim = 1:nsims,.packages = c('Rfast','SpatialExtremes')) %dopar% {
    set.seed(sim)
    output_local  <- veccia_app_local(Y = Y,locs=loc,ffnn_params = params_list_ensemble,
                                      K=15,iters=20000,burn=1000,thin=1,update=100,mean_delta = mean_deltas[sim])
} 
stopCluster(cl)

# Plots of the MCMC results
par(mfrow=c(2,4))
outHUC02 <- output_app[[1]]
hist(outHUC02[[2]][,1], main = expression(delta), xlab = '')
hist(outHUC02[[2]][,2], main = expression(rho), xlab = '')
hist(outHUC02[[2]][,3], main = expression(r), xlab = '')
hist(outHUC02[[2]][,4], main = expression(rho[2]), xlab = '')
outHUC02 <- output_app[[2]]
hist(outHUC02[[2]][,1], main = expression(delta), xlab = '')
hist(outHUC02[[2]][,2], main = expression(rho), xlab = '')
hist(outHUC02[[2]][,3], main = expression(r), xlab = '')
hist(outHUC02[[2]][,4], main = expression(rho[2]), xlab = '')

outHUC02 <- output_app[[1]]
plot(1:20000,outHUC02[[2]][,1], main = expression(delta), xlab = '',type = 'l')
plot(1:20000,outHUC02[[2]][,2], main = expression(rho), xlab = '',type = 'l')
plot(1:20000,outHUC02[[2]][,3], main = expression(r), xlab = '',type = 'l')
plot(1:20000,outHUC02[[2]][,4], main = expression(rho[2]), xlab = '',type = 'l')
outHUC02 <- output_app[[2]]
plot(1:20000,outHUC02[[2]][,1], main = expression(delta), xlab = '',type = 'l')
plot(1:20000,outHUC02[[2]][,2], main = expression(rho), xlab = '',type = 'l')
plot(1:20000,outHUC02[[2]][,3], main = expression(r), xlab = '',type = 'l')
plot(1:20000,outHUC02[[2]][,4], main = expression(rho[2]), xlab = '',type = 'l')
par(mfrow=c(1,1))