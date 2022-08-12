basis <- function(y,K,integral=FALSE){
    library(splines2)
    knots <- seq(1/(K-2), 1-1/(K-2), length=K-3)
    B     <- mSpline(y, knots = knots, 
                     Boundary.knots=c(0,1),
                     intercept=TRUE, degree = 2,
                     integral=integral)
    return(t(B))}

# Return softmax probabilities of an NN
probs_2l <- function(params,x){
    actfn    <- function(x){1/(1+exp(-x))}
    output  <- exp(params[[5]] %*% actfn(params[[3]] %*% actfn(params[[1]]%*%x + params[[2]]) + params[[4]]) + params[[6]])
    d       <- colsums(output)
    p       <- eachrow(output,d,'/')
    return(p)}

makeBA <- function(y0,X,theta,K,prev,id){
    BA      <- prev
    for(i in id){if(i>1){
        location = theta[i,1] + X*theta[i,2]
        Y        <- pgev(y0[[i]],location,exp(theta[i,3]),theta[i,4]) 
        BA[[i]]  <- basis(Y,K)
    }}
    return(BA)}

get_neighbor_theta = function(id,nn_id,theta,X_t){
    num_nn = length(na.omit(nn_id[id,]))
    neighs = nn_id[id,1:num_nn]
    theta_local = matrix(theta[neighs,],nrow = num_nn,ncol = 8)
    # print(theta_local)
    mu = matrix(nrow = num_nn,ncol = 50)
    sigma = matrix(exp(theta_local[,3]),nrow = num_nn,ncol = 50)
    xi = matrix(theta_local[,4],nrow = num_nn,ncol = 50)
    for(i in 1:num_nn){
        mu[i,] = theta_local[i,1] + theta_local[i,2]*X_t
    }
    return(list(mu,sigma,xi))
}

# Version where we extract parameters from NN and compute fitted values ourselves
makePR <- function(y1,X_t,params,theta,prev,id,nn_id){
    delta   <- pnorm(theta[1,5])
    if(delta<=0.5){pars <- params[[1]]}
    if(delta >0.5){pars <- params[[2]]}
    rho    <- 0.5*pnorm(theta[1,6])
    r    <- pnorm(theta[1,7])
    PR      <- prev
    for(i in id){if(i>1){
        # location = theta[i,1] + X_t*theta[i,2]
        gev_params = get_neighbor_theta(i,loc$nn_id,theta,X_t)
        Y       <- eva::pgev(y1[[i]],gev_params[[1]],gev_params[[2]],gev_params[[3]])
        X       <- features(delta=delta,rho=rho,r=r,Y=Y)
        PR[[i]] <- probs_2l(pars[[i]],X)
        # print(i)
    }}
    return(PR)}

log_like_app <- function(y0,X,theta,BA,PR,prev,id){
    ll    <- prev
    for(i in id){
        location = theta[i,1] + X*theta[i,2]
        ll[i] <- sum(dgev(y0[[i]],location,exp(theta[i,3]),theta[i,4],log=TRUE))+
            sum(log(colsums(BA[[i]]*PR[[i]])))
    }
    return(ll)}

inits <- function(y,X,show=FALSE){
    fit    <- ismev::gev.fit(xdat=as.vector(y),ydat=X,mul=1,siglink=exp,show=show)
    return(fit)}

dmvnorm <- function(y,mu,inv_cov,logdet){
    0.5*logdet - 0.5*t(y-mu)%*%inv_cov%*%(y-mu)   
}    

# mean_mu0=0;          sd_mu0=10;
# mean_mu1=0;          sd_mu1=10;
# mean_logsig=0;      sd_logsig=10;
# mean_xi = 0;        sd_xi = 0.25;
# mean_rho=0;        sd_rho=1;
# mean_delta=0;      sd_delta=1;
# mean_r=0;      sd_r=1;
# mean_rho2=-2;        sd_rho2=1;
# eps = 0.1;
# iters=5000;burn=1000;thin=1;update=1


veccia_app_local <- function(Y,locs,ffnn_params,K=15,
                             mean_mu0=0,          sd_mu0=10,
                             mean_mu1=0,          sd_mu1=10,
                             mean_logsig=0,      sd_logsig=10,
                             mean_xi = 0,        sd_xi = 0.25,
                             mean_rho=0,        sd_rho=1,
                             mean_delta=0,      sd_delta=1,
                             mean_r=0,      sd_r=1,
                             mean_rho2=-2,        sd_rho2=1,
                             eps = 0.1,
                             iters=5000,burn=1000,thin=1,update=1){
    
    tick <- proc.time()[3]
    
    # Extract QUINN features
    n      <- nrow(Y) # Number of sites
    m      <- ncol(Y) # Number of replications
    y0     <- list()
    y1     <- list()
    for(i in 1:n){
        y0[[i]] <- Y[i,]
        if(i>1){
            neighs  <- na.omit(locs$nn_id[i,])
            if(i==2){
                y1[[i]] <- matrix(Y[neighs,],nrow = 1)
            } else{
                y1[[i]] <- Y[neighs,]    
            }
            
        }
    }
    terms  <- list()
    for(i in 1:n){
        terms[[i]] <- i
        for(j in 1:n){if(i%in%locs$nn_id[j,]){
            terms[[i]] <- c(terms[[i]],j)
        }}
    }    
    
    # Bookkeeping 
    theta_mn <- c(mean_mu0, mean_mu1, mean_logsig, mean_xi, mean_delta, mean_rho, mean_r, mean_rho2)
    theta_sd <- c(  sd_mu0, mean_mu1, sd_logsig,   sd_xi,   sd_delta,   sd_rho,   sd_r,   sd_rho2)
    p        <- length(theta_mn)
    
    years   <- 1972:2021
    X_t     <- 0.1*(years-mean(years))
    theta   <- matrix(c(0,0,0,0,mean_delta,mean_rho,mean_r,mean_rho2),nrow = n,ncol = p,byrow = T)
    for(j in 1:n){
        init_fit = inits(y=Y[j,],matrix(X_t,50,1))
        theta[j,1:4] = init_fit$mle
    }
    
    mu          <- theta[,1:4]
    # for(j in 2:n){
    #     print(j)
    #     OK <- max(pgev(y0[[j]],theta[j,1] + X_t*theta[j,2],exp(theta[j,3]),theta[j,4]))<1 & min(pgev(y1[[j]],theta[j,1] + X_t*theta[j,2],exp(theta[j,3]),theta[j,4]))>0
    #     while(!OK){
    #         print('update')
    #         theta[j,3] <- theta[j,3] + 0.05
    #         OK         <- max(pgev(y0[[j]],theta[j,1] + X_t*theta[j,2],exp(theta[j,3]),theta[j,4]))<1 & min(pgev(y0[[j]],theta[j,1] + X_t*theta[j,2],exp(theta[j,3]),theta[j,4]))>0
    #     }
    # }
    # theta[18,3] = log(2*exp(theta[18,3]))
    #Initialize the SVC parameters
    tau         <- rep(1,4)
    b1          <- rep(0,4)
    sig         <- rep(0.1,4)
    rho2        <- mean_rho2
    DIST        <- fields::rdist(locs$locs)
    SIGMA       <- exp(-DIST/exp(rho2))
    OMEGA       <- solve(SIGMA)
    logdet      <- determinant(OMEGA)$modulus
    ones        <- matrix(1,nrow = n, ncol = 1)

    


    # Initialize remaining variables
    params   <- ffnn_params
    BA       <- list()
    PR       <- list()
    BA[[1]]  <- diag(2)
    PR[[1]]  <- diag(2)
    BA       <- makeBA(y0,X_t,theta,K,prev=BA,id=1:n)
    PR       <- makePR(y1,X_t,params,theta,prev=PR,id=1:n,nn_id = loc$nn_id)
    sapply(PR,function(x)sum(is.nan(x)))
    curll    <- log_like_app(y0,X_t,theta,BA,PR,prev=rep(0,n),id=1:n)
    curll

    # Keep track of stuff
    keep_theta_gev          <- array(0,c(iters,n,4))
    keep_theta_sp           <- matrix(0,iters,4)
    colnames(keep_theta_sp) <- c("delta","rho","r","rho2")
    att  <- acc <- MH <- c(1,1,1,0.1,.1,.1,.1,.1)
    P_MH <- t(chol(cor(mu)))
    
    # GO!!!
    
    for(iter in 1:iters){
        for(ttt in 1:thin){
            # print(iter)
            # aa = proc.time()[3]
            ################################:
            ####   UPDATE PARAMETERS   #####:
            ################################:
            # Updated GEV parameters
            for(i in 1:n){
                ti       <- terms[[i]] 
                att[1]   <- att[1] + 1
                cant     <- theta
                cant[i,1:4] <- theta[i,1:4] + MH[1]*tau[1:4]*as.vector(P_MH%*%rnorm(4))
                canBA    <- makeBA(y0,X_t,cant,K,prev=BA,id=ti)
                canPR    <- makePR(y1,X_t,params,cant,prev=PR,id=ti)
                canll    <- log_like_app(y0,X=X_t,theta=cant,BA=canBA,PR=canPR,prev=curll,id=ti)
                mhprob   <- sum(canll[ti]-curll[ti])+
                    sum(dnorm(cant[i,1:4],mu[i,1:4],tau[1:4],log=TRUE))-
                    sum(dnorm(theta[i,1:4],mu[i,1:4],tau[1:4],log=TRUE))
                if(!is.na(mhprob)){if(log(runif(1))<mhprob){
                    acc[1]    <- acc[1] + 1
                    theta[i,] <- cant[i,]
                    BA  <- canBA
                    PR  <- canPR
                    curll[ti] <- canll[ti]
                }}  
            }
            # print(paste(sum(is.nan(canll))),sum(is.nan(curll)))
            # print(mhprob)
            # Update random effects distribution parameters
    
            for(j in 1:4){
                R       <- theta[,j] - mu[,j]
                tau[j]  <- 1/sqrt(rgamma(1,0.5*n+eps,0.5*sum(R^2)+eps))

                R1      <- mu[,j] - b1[j]
                R2      <- t(R1) %*% OMEGA %*% R1
                sig[j]  <- 1/sqrt(rgamma(1,0.5*n+eps, 0.5*R2+eps))

                VINV    <- (1/100) + sum(OMEGA)/(sig[j]^2)
                MMM     <- (1/sig[j]^2)*t(ones) %*% OMEGA %*% mu[,j]
                b1[j]   <- rnorm(1,MMM/VINV, 1/sqrt(VINV))

                VINV2   <- solve((1/tau[j]^2)*diag(n) + (1/sig[j]^2)*OMEGA)
                MMM2    <- (1/tau[j]^2)*theta[,j] + (b1[j]/sig[j]^2)*(OMEGA %*% ones)
                # mu[,j]  <- MASS::mvrnorm(1,VINV2%*%MMM2, VINV2)
                mu[,j]  <- mvnfast::rmvn(1,VINV2%*%MMM2, VINV2)

            }
            
#Update log range
            att[p]    <- att[p] + 1
            canrho2   <- rnorm(1,rho2, MH[p])
            
            canSIGMA  <- exp(-DIST/exp(canrho2))
            canOMEGA  <- solve(canSIGMA)
            canlogdet <- determinant(canOMEGA)$modulus            
            MHprob    <- dnorm(canrho2,mean_rho2,sd_rho2,log=TRUE)-
                         dnorm(   rho2,mean_rho2,sd_rho2,log=TRUE)
            for(j in 1:4){
                MHprob <- MHprob + dmvnorm(mu[,j],b1[j],canOMEGA/(sig[j]^2),canlogdet)-
                                   dmvnorm(mu[,j],b1[j],OMEGA/(sig[j]^2),logdet)
            }
            # print(paste(canrho2,MHprob))
            if(log(runif(1))<MHprob){
                # print('success')
                acc[p]      <- acc[p] + 1
                theta[,8] <- rho2        <- canrho2
                SIGMA       <- canSIGMA
                OMEGA       <- canOMEGA
                logdet      <- canlogdet
            }
            
            # Updates spatial dependence parameters
            
            for(j in 5:(p-1)){
                att[j]   <- att[j] + 1
                cant     <- theta
                cant[,j] <- rnorm(1,theta[1,j],MH[j])
                canPR    <- makePR(y1,X_t,params,cant,prev=PR,id=1:n)
                canll   <- log_like_app(y0,X_t,cant,BA,canPR,prev=curll,id=1:n)
                mhprob  <- sum(canll-curll)+
                    dnorm( cant[1,j],theta_mn[j],theta_sd[j],log=TRUE)-
                    dnorm(theta[1,j],theta_mn[j],theta_sd[j],log=TRUE)
                # print(paste(sum(is.nan(canll))),sum(is.nan(curll)))
                # print(mhprob)
                if(!is.na(mhprob)){if(log(runif(1))<mhprob){
                    acc[j] <- acc[j] + 1
                    theta  <- cant
                    PR     <- canPR
                    curll  <- canll
                }}  
            }
            
            # Tuning
            for(j in 1:length(att)){if(att[j]>50 & iter<burn){
                if(acc[j]/att[j]<0.2){MH[j] <- MH[j]*0.8}
                if(acc[j]/att[j]>0.6){MH[j] <- MH[j]*1.2}
                acc[j] <- att[j] <- 0
            }}
            
        } # end thin
        
        ##############################################:
        #####        KEEP TRACK OF STUFF       #######:
        ##############################################:
        
        keep_theta_gev[iter,,] <- theta[,1:4]
        keep_theta_gev[iter,,3] <- exp(theta[,3])
        keep_theta_sp[iter,]   <- c(pnorm(theta[1,5]),0.5*pnorm(theta[1,6]),pnorm(theta[1,7]),exp(theta[1,8]))
        
        ##############################################:
        #####       PLOT RESULTS SO FAR        #######:
        ##############################################:
        if(update>0 & iter%%update==0){
            print(paste(iter,pnorm(theta[1,5]),0.5*pnorm(theta[1,6]),pnorm(theta[1,7]),exp(theta[1,8])))
            par(mfrow=c(2,4))
            for(j in 1:4){
                plot(keep_theta_gev[1:iter,21,j],type="l",
                     xlab="MCMC iteration",ylab="Sample")
            }
            for(j in 1:4){
                plot(keep_theta_sp[1:iter,j],type="l",
                     xlab="MCMC iteration",ylab="Sample",
                     main=colnames(keep_theta_sp)[j])
            }
        }
        # print(proc.time()[3]-aa)
    }   
    tock <- proc.time()[3]
    out  <- list(theta_gev=keep_theta_gev,theta_sp=keep_theta_sp,time=tock-tick)
    return(out)}
