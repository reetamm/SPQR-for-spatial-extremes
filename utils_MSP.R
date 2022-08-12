library(pracma)
# library(Rmpfr)
make.coordinates = function(h){
  coord = h/sqrt(2)
  a = matrix(c(0,0,coord,coord),nrow = 2,byrow = T)
  return(a)
}

br_to_exp = function(x){
  a = -log(1-exp(-1/x))
  return(a)
}

gp_to_exp = function(x,nugget=1){
  a = -log(1-pnorm(x,sd=sqrt(nugget+1)))
  return(a)
}

br_to_unif = function(x){
  a = exp(-1/x)
  return(a)
}

chi_h = function(data,u){
  F_j = data[,1]
  F_k = data[,2]
  freq_joint = sum(ifelse(F_j>u & F_k > u,1,0))
  freq_margin = sum(ifelse(F_k>u,1,0))
  return(freq_joint/freq_margin)
}

phypoexp2 = function(data,rate1,rate2){
  cdf = 1 - (rate2/(rate2-rate1))*exp(-rate1*data) + (rate1/(rate2-rate1))*exp(-rate2*data)
  return(cdf)
}

phypoexp2_log = function(data,rate1,rate2){
    cdf = 1 - exp(log(rate2/(rate2-rate1))-rate1*data) + exp(log(rate1/(rate2-rate1))-rate2*data)
    return(cdf)
}

extremal_coeff = function(h, range, smooth){
    semivar = 0.5*(h/range)^smooth
    vartheta = 2*pnorm(sqrt(semivar))
    return(vartheta)
}

msp_cdf = function(h, range, smooth,x_1,x_2){
  semivar = 0.5*(h/range)^smooth
  a_12 = sqrt(2*semivar)
  V_12 = (1/x_1)*pnorm(0.5*a_12 - log(x_1/x_2)/a_12) + 
    (1/x_2)*pnorm(0.5*a_12 - log(x_2/x_1)/a_12)
  return(exp(-V_12))
}


chi_h_common = function(delta, tail){
    exp1 = exp(-tail*(1-2*delta)/(delta*(1-delta)))
    exp2 = exp(-tail/(1-delta))
    numerator = ((1-delta)/(3*delta-1))*(exp1-exp2) + exp1
    denominator = ((1-delta)/(1-2*delta)) - (delta/(1-2*delta))*exp1
    return(numerator/denominator)
}

g_r_inv = function(x){
    x1 = log(1-exp(-x))
    x2 = -1/x1
    return(x2)
}


testfn = function(h, range, smooth, delta, tail){
    upper.limit = tail/(1-delta)
    f = function(x_1,x_2,h, range, smooth,delta, tail){
        r_1 = (tail - (1-delta)*x_1)/delta
        r_2 = (tail - (1-delta)*x_2)/delta
        # print(paste(g_r_inv(r_1),g_r_inv(r_2)))
        # print(paste(r_1,r_2))
        arg = log(msp_cdf(h, range, smooth, g_r_inv(r_1),g_r_inv(r_2))) -x_1 - x_2
        return(exp(arg))
    }
    # f2 = function(x,h, range, smooth,delta, tail){
    #     r_1 = (tail - (1-delta)*x[1])/delta
    #     r_2 = (tail - (1-delta)*x[2])/delta
    #     arg = msp_cdf(h, range, smooth, g_r_inv(r_1),g_r_inv(r_2))*exp(-x[1] - x[2])
    #     return(arg)
    # }
    t5 = integral2(f,0,upper.limit,0,upper.limit,singular=F,h=h,range=range,smooth=smooth,delta = delta, tail=tail)
    # t5 = pcubature(f2,c(0,0),c(upper.limit,upper.limit),h=h,range=range,smooth=smooth,delta = delta, tail=tail)
    t1 = exp(-tail/delta)*(delta/(1-2*delta))*(exp(tail*(1-2*delta)/(delta*(1-delta))) - 1)
    t2 = t1
    t3 = 1 + mpfr(exp(-2*tail/(1-delta)) - 2*exp(-tail/(1-delta)),300)
    # t3 = 1 + exp(-2*tail/(1-delta)) - 2*exp(-tail/(1-delta))
    t4 = exp(-2*tail/(1-delta))
    # if(t5$Q==1){
    #     print(1)
    # } else{
    #     print(log(t5$Q))
    # }
    # if(log(t5$Q) > -5e-8)
    #     t5$Q = 1
    numerator = t1+t2-t3+t4+t5$Q
    denominator = ((1-delta)/(1-2*delta))*exp(-tail/(1-delta)) - (delta/(1-2*delta))*exp(-tail/delta)
    return(numerator/denominator)
    # return(c(t1,t2,t3,t4,t5$Q,denominator))
    # return(c(numerator,denominator,t5-t3))
    # return(t5)
}
expit<-function(x){1/(1+exp(-x))}