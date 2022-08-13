rm(list = ls())
library(fields)
library(torch)
library(SPQR)
library(SpatialExtremes)
library(doParallel)
source("gen_data.R",local = T)
source("utils_MSP.R",local = T)
load("HCDN/HCDN_annual_max.RData")

# Switches to GPU if CUDA is set up.
# This is a flag in the SPQR function
if(cuda_is_available())
    use.gpu = T
if(!cuda_is_available())
    use.gpu = F

trunc.data = data_prep(Y,s)
Y = trunc.data$Y
s = trunc.data$s

n       <- nrow(s) # Number of locations
loc     <- create_locations(NS = 15, m = n, locs = s) # NS is number of neighbors

n_train     <- 200000 # Training data size for each local SPQR
lr          <- 0.01 # Learning Rate
epochs      <- 50 # Number of training epochs
n.knots     <- 15 # Number of output knots
n.hidden    <- c(30,20) # Number of neurons in each hidden layer

features    <- function(delta,rho,r,Y){qnorm(rbind(delta,2*rho,r,Y))} # Function to create covariate matrix

# When creating a cluster, make sure that the first argument is less than the number
# of cores in your system. Unless you are sure, it is recommended to compute
# sequentially. Uncomment registerDoParallel and comment out registerDoSEQ only if you are sure your 
# computer can handle it

cl          <- makeCluster(10,outfile = 'Log.txt')
                                                    
# registerDoParallel(cl)
registerDoSEQ()

fits <- foreach(id=2:489, .packages = c('SPQR','torch','SpatialExtremes')) %dopar% {
    mle.control     <- list(batch.size = 1000, epochs = epochs, lr = lr,
                            use.GPU = use.gpu, early.stopping.epochs = 20, save.name=paste0("mspgp_fit.pt"))
    tic <- proc.time()
    min_delta <- 0.01
    max_delta <- 0.99
    set.seed(919*id + 1)
    nn          <- c(id,na.omit(loc$nn_id[id,])) # Each local SPQR will generate data for only the location and its conditioning set
    dat         <- gen_data_evp_local(n = n_train, s = as.matrix(loc$locs[nn,]),
                                     min_delta = min_delta, max_delta = max_delta)
    rho         <- dat$rho
    delta       <- dat$delta
    r           <- dat$r
    Y           <- dat$Y
    modelname   <- paste('ffnn', id, '1', sep = '_')
    mle.control$save.name <- paste0("mspgp_fit_st", id, "batch 1.pt")
    y           <- Y[1,]
    X           <- t(features(delta = delta, rho = rho, r = r, Y = Y[-1,]))   
    fit         <- SPQR(X = X, Y = y, n.knots = n.knots, n.hidden = n.hidden, activation = "sigmoid",
                      method = "MLE", control = mle.control, normalize = FALSE, verbose = F)
    save.SPQR(fit,modelname, path = 'EVP_HUC02')
    toc <- proc.time()
    timetaken <- (toc-tic)[3]
    print(paste('location', id, 'batch 1 fitted in', timetaken, 'seconds'))
}
stopCluster(cl)