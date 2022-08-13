logit <- function(x){log(x/(1-x))}

# Prep the HUC02 data for analysis
data_prep = function(Y,s){
    Y       <- Y[,23:72] # Select last 50 years
    year    <- year[23:72] # Select last 50 years
    nmiss   <- rowSums(is.na(Y)) # Missingness at stations
    Y       <- Y[nmiss==0,] # Select stations with complete data
    s       <- s[nmiss==0,] # Coordinates of stations with complete data
    
    # Scale the domain to (0,1)
    s[,1]   <- s[,1]-min(s[,1])
    s[,2]   <- s[,2]-min(s[,2])
    s_scale <- max(s)
    s       <- s/s_scale
    
    Y[Y<1]  <- 1 # Set minimum value to 1
    Y       <- log(Y) # log transform
    out = list(Y=Y,s=s,nmiss=nmiss,s_scale = s_scale)
    return(out)
}

# Sample spatial locations and compute neighbors
create_locations_evp <- function(k=10,n=50,seed=1,min_dist=0.025){
    set.seed(seed)
    locs     <- cbind(runif(n),runif(n))
    mind     <- min(dist(locs))
    while(mind < min_dist){
        d        <- as.matrix(dist(locs))
        diag(d)  <- Inf
        d        <- which.min(apply(d,1,min))[1]
        locs[d,] <- runif(2)
        mind     <- min(dist(locs))
    }
    locs       <- locs[order(rowSums(locs)),]
    nn_dist    <- array(Inf,c(k,2,n))
    nn_id      <- matrix(NA,n,k)
    dd         <- as.matrix(dist(locs))      
    for(i in 2:n){
        if(i<=(k+1)){
            nn_id[i,1:(i-1)]     <- 1:(i-1)
            nn_dist[1:(i-1),1,i] <- locs[1:(i-1),1]-locs[i,1]
            nn_dist[1:(i-1),2,i] <- locs[1:(i-1),2]-locs[i,2]
        }
        if(i>(k+1)){
            dtrunc        <- ifelse(1:n <i,dd[i,],Inf)
            nn_id[i,]     <- which(rank(dtrunc)<=k)
            nn_dist[,1,i] <- locs[nn_id[i,],1]-locs[i,1]
            nn_dist[,2,i] <- locs[nn_id[i,],2]-locs[i,2]
        }
    }
    nn_dist[nn_dist==Inf] <- 10 
    for(i in 2:n){
        dist <- nn_dist[,1,i]^2 + nn_dist[,2,i]^2
        nn_order <- order(dist)
        nn_dist[,1,i] <- nn_dist[nn_order,1,i]
        nn_dist[,2,i] <- nn_dist[nn_order,2,i]
        nn_id[i,] <- nn_id[i,nn_order]
    }
    out     <- list(locs=locs,nn_dist=nn_dist,nn_id=nn_id,k=k)
    return(out)
    }


create_locations <- function(NS=10,m=100,seed=1,locs=NULL){
  rs <- .Random.seed 
  set.seed(seed)
  if(is.null(locs))
      locs  <- cbind(runif(m),runif(m))
  locs      <- locs[order(rowSums(locs)),]
  nn_dist   <- array(Inf,c(NS,2,m))
  nn_id     <- matrix(NA,m,NS)
  dd        <- as.matrix(dist(locs))      
  for(i in 2:m){
   if(i<=(NS+1)){
      nn_id[i,1:(i-1)]     <- 1:(i-1)
      nn_dist[1:(i-1),1,i] <- locs[1:(i-1),1]-locs[i,1]
      nn_dist[1:(i-1),2,i] <- locs[1:(i-1),2]-locs[i,2]
   }
    if(i>(NS+1)){
     dtrunc        <- ifelse(1:m <i,dd[i,],Inf)
     nn_id[i,]     <- which(rank(dtrunc)<=NS)
     nn_dist[,1,i] <- locs[nn_id[i,],1]-locs[i,1]
     nn_dist[,2,i] <- locs[nn_id[i,],2]-locs[i,2]
    }
  }
  nn_dist[nn_dist==Inf] <- 10 
  for(i in 2:m){
      dist <- nn_dist[,1,i]^2 + nn_dist[,2,i]^2
      nn_order <- order(dist)
      nn_dist[,1,i] <- nn_dist[nn_order,1,i]
      nn_dist[,2,i] <- nn_dist[nn_order,2,i]
      nn_id[i,]     <- nn_id[i,nn_order]
  }
  out     <- list(locs=locs,nn_dist=nn_dist,nn_id=nn_id,k=NS)
  .Random.seed <-rs
  return(out)
  }

rhoW2rhoR <- function(rhoW,alpha){
    rhoW*(log(20)/(4*(qnorm(1-0.05/2)^2)))^(1/alpha)  
}

# rhoR and rhoW related deterministically, also presence of nugget
gen_data_evp_local <- function(n,s,min_rho=0.01,max_rho=0.50,min_delta=0.01,max_delta=0.99,min_r = 0.01, max_r = 0.99){
    rho     <- runif(n,min_rho,max_rho)
    delta   <- runif(n,min_delta,max_delta)
    r       <- runif(n, min_r, max_r)
    Y       <- matrix(nrow = nrow(s),ncol = n)
    for(i in 1:n){
        br_local <- t(rmaxstab(1, s, "brown", range = rhoW2rhoR(rho[i],1), smooth = 1))
        gev_mar <- matrix(rgev(nrow(s), 1, 1, 1), ncol = 1)
        tmp <- cbind(r[i]*br_local,(1-r[i])*gev_mar)
        R_s <- apply(tmp,1,max)
        gp_local<- t(rgp(1, s, "powexp", nugget = r[i], range = rho[i], smooth = 1))
        br_exp  <- br_to_exp(R_s)
        gp_exp  <- gp_to_exp(gp_local,nugget = r[i])
        y       <- delta[i]*br_exp + (1-delta[i])*gp_exp
        Y[,i]   <- phypoexp2(y, rate1 = 1/delta[i], rate2 = 1/(1-delta[i]))
    }
    out <- list(Y=Y,rho=rho,delta=delta,r=r)
    return(out)
    }


# Construct the feature for deep learning
features.evp.local <- function(delta,rhoR,rhoW,Y){
    x     <- rbind(delta,rhoR,rhoW,Y)
    return(x)
    }

# Extract FFNN weights from an SPQR object
get.nn.params <- function(fitted.obj){
    a <- fitted.obj$model$parameters
    ffnn_params <- list()
    for(j in 1:length(a)){
        ffnn_params[[j]] <- torch::as_array(a[[j]])
    }
    return(ffnn_params)
}

# Misc functions for plotting
bp <- function(x1,truth,main="",fx=function(x){x}){
    boxplot(fx(x1),outline=FALSE,names=c("Approx"),main=main)
    abline(fx(truth),0,col=2)
}

fx=function(x){x}

hist_mean <- function(x1,truth,main="",fx=function(x){x}){
    hist(fx(x1),main=main)
    abline(v=fx(truth),col=2)
}

###########################################
# Functions used for the simulation studies
###########################################

# Kriging equations
krige <- function(s,Y,r,rho){
    s0   <- rbind(0,s)
    d    <- as.matrix(dist(s0))
    d0   <- d[1,-1]
    d    <- d[-1,-1]
    P22  <- solve((1-r)*diag(nrow(d)) + r*exp(-d/rho))
    S12  <- r*exp(-d0/rho)
    PS   <- t(S12)%*%P22
    M    <- sum(PS*Y)
    V    <- 1 - sum(PS*S12)
    out  <- list(mean=M,sd=sqrt(V))
    return(out)
    }

# Gaussian conditional PDF
dYgivenX <- function(y,r,rho,Y,s){
    k    <- krige(s,Y,r,rho)
    d    <- dnorm(qnorm(y),k$mean,k$sd)/dnorm(qnorm(y))
    return(d)
    }

# Gaussian conditional quantile function
qYgivenX <- function(tau,r,rho,Y,s){
    k    <- krige(s,Y,r,rho)
    q    <- pnorm(qnorm(tau,k$mean,k$sd))
    return(q)
    }

# Sample training data
gen_data_gp_global <- function(n,ns,loc,ar=1,br=1,min_rho=0.01,max_rho=2){
    m      <- dim(loc$locs)[1]
    r      <- rbeta(n,ar,br)
    rho    <- runif(n,min_rho,max_rho)
    id     <- sample(2:m,n,replace=TRUE)
    s      <- loc$nn_dist[,,id]
    Y      <- matrix(0,ns+1,n)
    for(i in 1:n){
        ss         <- rbind(0,s[,,i])
        d          <- as.matrix(dist(ss))
        COV        <- (1-r[i])*diag(nrow(d)) + r[i]*exp(-d/rho[i]) 
        Y[,i]      <- t(chol(COV))%*%rnorm(ns+1) 
    }
    out <- list(Y=Y,s=s,r=r,rho=rho,id=id)
    return(out)
    }

gen_data_gp_local <- function(n,s,ar=1,br=1,min_rho=0.01,max_rho=2){
    m      <- nrow(s) 
    d      <- fields::rdist(s)
    r      <- rbeta(n,ar,br)
    rho    <- runif(n,min_rho,max_rho)
    Y      <- matrix(0,m,n)
    for(i in 1:n){
        COV        <- (1-r[i])*diag(m) + r[i]*exp(-d/rho[i])
        Y[,i]      <- t(chol(COV))%*%rnorm(m)
    }
    out <- list(Y=Y,r=r,rho=rho)
    return(out)
    }

# Independent rhoR and rhoW true values
gen_data_evp_local1 <- function(n,s,min_alpha=1,max_alpha=2,min_rho=0.01,max_rho=1,min_delta=0,max_delta=1,verbose=F){
    m       <- nrow(s) 
    alpha   <- 1+0*runif(n,min_alpha,max_alpha)
    rhoR    <- runif(n,min_rho,max_rho)
    rhoW    <- runif(n,min_rho,max_rho)
    delta   <- runif(n,min_delta,max_delta)
    Y       <- matrix(0,m,n)
    for(i in 1:n){
        br_local <- rmaxstab(1, s, "brown", nugget = 0, range = rhoR[i], smooth = alpha[i])
        gp_local <- rgp(1, s, "powexp", nugget = 0, range = rhoW[i], smooth = alpha[i])
        br_exp   <- br_to_exp(br_local)
        gp_exp   <- gp_to_exp(gp_local)
        Y[,i]    <- delta[i]*br_exp + (1-delta[i])*gp_exp
        Y[,i]    <- phypoexp2(Y[,i], rate1 = 1/delta[i], rate2 = 1/(1-delta[i]))
        if(verbose==T)
            print(paste(i,'of',n))
    }
    out <- list(Y=Y, rhoR=rhoR, rhoW=rhoW, delta=delta)
    return(out)
    }

# Shared value of rhoR and rhoW
gen_data_evp_local2 <- function(n,s,min_rho=0.01,max_rho=0.99,min_delta=0.01,max_delta=0.99){
    m      <- nrow(s) 
    rho    <- runif(n,min_rho,max_rho)
    delta  <- runif(n,min_delta,max_delta)
    Y      <- matrix(0,m,n)
    for(i in 1:n){
        br_local <- rmaxstab(1, s, "brown", nugget = 0, range = rho[i], smooth = 1)
        gp_local <- rgp(1, s, "powexp", nugget = 0, range = rho[i], smooth = 1)
        br_exp   <- br_to_exp(br_local)
        gp_exp   <- gp_to_exp(gp_local)
        Y[,i]    <- delta[i]*br_exp + (1-delta[i])*gp_exp
        Y[,i]    <- phypoexp2(Y[,i], rate1 = 1/delta[i], rate2 = 1/(1-delta[i]))
    }
    out <- list(Y=Y,rho=rho,delta=delta)
    return(out)
    }

# Construct the feature for deep learning
features.gp.global <- function(s,Y,rho,r){
    x <- rbind(log(r/(1-r)),Y,
               sweep(s[,1,],2,rho,"/"),
               sweep(s[,2,],2,rho,"/"))
    return(x)
    }

features.gp.local <- function(Y,rho,r){
    x <- rbind(log(r/(1-r)),log(rho),Y)
    return(x)
    }