library(SpatialExtremes)
library(ggplot2)
library(graphics)
library(foreach)
library(doParallel)

cores = detectCores()
cores
cl <- makeCluster(cores[1]-1) # not to overload your computer
registerDoParallel(cl)
# registerDoSEQ()
distances = seq(0.1,2,by=0.1)
source("./utils_MSP.R")
source("./gen_data.R")

coord = make.coordinates(2)
rho = 0.5
rhoR = rhoW2rhoR(rho,alpha = 1)
n=10000000
br_list_HW = foreach(h = distances, .packages = 'SpatialExtremes') %dopar% {
  coord = make.coordinates(h)
  # fit                <- veccia_GP(Y,locs,burn=burn,iters=iters,update=update)
  print(h)
  br_local <- rmaxstab(n, coord, "brown", nugget = 0, range = rhoR, smooth = 1)
} 

gp_list_HW = foreach(h = distances, .packages = 'SpatialExtremes') %dopar% {
  coord = make.coordinates(h)
  # fit                <- veccia_GP(Y,locs,burn=burn,iters=iters,update=update)
  print(h)
  gp_local <- rgp(n, coord, "powexp", nugget = 0, range = rho, smooth = 1)
} 
(-log(.4))^2
distances[8]

br1_data = br_list_HW[[8]]
gp1_data = gp_list_HW[[8]]
# br1_data[,2] = br1_data[,1]
br1_exp = br_to_exp(br1_data)
hist(br1_exp[,1])
br1_exp[,2] = br1_exp[,1]
br1_unif = br_to_unif(br1_data)
hist(br1_unif[,1])

gp1_exp = gp_to_exp(gp1_data)
hist(gp1_exp[,1])

u.sequence = seq(from = 0.900001, to = 0.999999, length = 10000)
chi_h_values = matrix(NA,10000,4)
delta.vector = c(.2,.4,.6,.8)

for(i in 1:4){
  delta = delta.vector[i]
  print(delta)
  mix.data = delta*br1_exp + (1-delta)*gp1_exp
  mix.data.cdf = phypoexp2(mix.data,rate1 = 1/delta,rate2 = 1/(1-delta))
  keep.locs = which(mix.data.cdf[,1]>.8 & mix.data.cdf[,2]>.8)
  mix.data.cdf = mix.data.cdf[keep.locs,]
  if(max(mix.data.cdf[,1])>max(mix.data.cdf[,2])){
      tmp = cbind(mix.data.cdf[,2],mix.data.cdf[,1])
      mix.data.cdf = tmp
  }
  for(j in 1:length(u.sequence)){
    u = u.sequence[j]
    chi_h_values[j,i] = chi_h(mix.data.cdf,u)
    print(paste(i,j,chi_h_values[j,i]))
  }
}

par(mfrow = c(1,1))
locs = which(u.sequence>.9)
a=matplot(u.sequence[locs],chi_h_values[locs,],'l',lty = 1, lwd = 2) +title('empirical common R')
legend(.95, .6, c("d = 0.2", "d = 0.4",
               "d = 0.8"),
       pch = "123",col=1:3)
tail(chi_h_values)
chi_h_values = chi_h_values[1:9999,]
locs = locs[1:9999]
all_chi = c(chi_h_values[locs,1],chi_h_values[locs,2],chi_h_values[locs,3],chi_h_values[locs,4])
howmany = length(all_chi)/4
deltas = factor(c(0.2,0.4,.6,0.8))
chi_df = data.frame(chi = all_chi, u = rep(u.sequence[locs],4),delta = rep(deltas,each=howmany) )
ggplot(chi_df,aes(x=u,y=chi,col=delta)) + geom_line(size=1) + theme_bw() +
    labs(x='Quantile (u)',y=expression(chi[u](h)),colour=expression(delta)) + ylim(0,1) +
    theme(legend.position = 'bottom',
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
          legend.spacing.x = unit(.7, 'cm'))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )+
    scale_color_manual(values = c('Black','Red','Green','Blue')) + 
    geom_vline(xintercept = 1, linetype='dotted',color='black',size = 1)



# y = sort(mix.data[,1])
# y_cdf = sort(mix.data.cdf[,1])
# tmp1 = chi_h_common(.1,tail=y)
# tmp2 = chi_h_common(.4,tail=y)
# tmp3 = chi_h_common(.7,tail=y)
# tmp4 = chi_h_common(.9,tail=y)
# locs = which(y_cdf>.8)
# chi_common = cbind(tmp1,tmp2,tmp3,tmp4)
# b=matplot(y_cdf[locs],chi_common[locs,],'l',lty = 1, lwd = 2) +title('theoretical common R')
# legend(.95, .6, c("d = 0.1", "d = 0.4",
#                   "d = 0.7", "d = 0.9"),
#        pch = "1234",col=1:4)
# locs2 = which(u.sequence<.99)
# d=matplot(y_cdf[locs2],chi_common[locs2,]-chi_h_values[locs2,],'l',lty = 1, lwd = 2) +title('theoretical - empirical')
# legend(0, -0.05, c("d = 0.1", "d = 0.4",
#                   "d = 0.7", "d = 0.9"),
#        pch = "1234",col=1:4)


chi_h_values2 = matrix(NA,100,12)
tmp = matrix(NA,100,4)
chi_h_values2 = cbind(chi_h_values2,tmp)
# a = c(1,1.5,1.95)
a=1
distances = seq(0.01,2,length=100)
for(k in 1:1){
  br_list2 = foreach(h = distances, .packages = 'SpatialExtremes') %dopar% {
    coord = make.coordinates(h)
    # fit                <- veccia_GP(Y,locs,burn=burn,iters=iters,update=update)
    br_local <- rmaxstab(1000000, coord, "brown", nugget = 0, range = rhoR, smooth = a[k])
  }

  gp_list2 = foreach(h = distances, .packages = 'SpatialExtremes') %dopar% {
    coord = make.coordinates(h)
    # fit                <- veccia_GP(Y,locs,burn=burn,iters=iters,update=update)
    gp_local <- rgp(1000000, coord, "powexp", nugget = 0, range = rho, smooth = a[k])
  }

  print(paste('Generated dataset with alpha =',a[k]))
for(i in 1:100){
  print(i)
  br_exp = br_to_exp(br_list2[[i]])
  gp_exp = gp_to_exp(gp_list2[[i]])
  for(j in 1:4){
    delta = delta.vector[j]
    mix.data = delta*br_exp + (1-delta)*gp_exp
    mix.data.cdf = phypoexp2(mix.data,rate1 = 1/delta,rate2 = 1/(1-delta))
    column = (k-1)*4 + j
    chi_h_values2[i,column] = chi_h(mix.data.cdf,0.99)
  }
}
}
par(mfrow = c(1,1))


matplot(distances,chi_h_values2[,1:4],'l',lty = 1,lwd=2) + title('alpha = 1.0')
legend(6, .8, c("d = 0.1", "d = 0.4",
                "d = 0.7", "d = 0.9"),
       pch = "12345",col=1:4)

# matplot(distances,chi_h_values2[,5:8],'l',lty = 1,lwd=2)+ title('alpha = 1.5')
# legend(6, .8, c("d = 0.1", "d = 0.4",
#                 "d = 0.7", "d = 0.9"),
#        pch = "12345",col=1:4)
# 
# matplot(distances,chi_h_values2[,9:12],'l',lty = 1,lwd=2)+ title('alpha = 1.95')
# legend(6, .8, c("d = 0.1", "d = 0.4",
#                 "d = 0.7", "d = 0.9"),
#        pch = "12345",col=1:4)
# 

tmp1 = matrix(chi_h_values2[,1:4])
tmp2 = matrix(chi_h_values2[,5:8])
tmp3 = matrix(chi_h_values2[,9:12])

chi_df_2 = data.frame(chi1 = tmp1,chi2 = tmp2,chi3 = tmp3,distances=rep(distances,4),delta = rep(deltas,each=100))
ggplot(chi_df_2,aes(x=distances,y=chi1,col=delta)) + geom_line(size=1) +#geom_smooth(size=1.5) + 
    theme_bw() +
    labs(x='Distance (h)',y=expression(chi[u](h)),colour=expression(delta)) + ylim(0,.8) +
    theme(legend.position = 'bottom',
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
          legend.spacing.x = unit(.7, 'cm'))+
    guides(color = guide_legend(override.aes = list(size = 2) ) ) +
    scale_color_manual(values = c('Black','Red','Green','Blue')) 

ggplot(chi_df_2,aes(x=distances,y=chi2,col=delta)) + geom_line(size=1) + #geom_smooth(size=1.5) + 
    theme_bw() +
    labs(x='Distance (h)',y=expression(chi[u](h)),colour=expression(delta)) + ylim(0,.8) +
    theme(legend.position = 'bottom',
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
          legend.spacing.x = unit(.7, 'cm'))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )+
    scale_color_manual(values = c('Black','Red','Green','Blue'))+
    annotate('text',x=8.75,y=0.73,size=8,label = expression(paste(alpha,' = 1.5')))

ggplot(chi_df_2,aes(x=distances,y=chi3,col=delta)) + geom_line(size=1) +#geom_smooth(size=1.5) + 
    theme_bw() +
    labs(x='Distance (h)',y=expression(chi[u](h)),colour=expression(delta)) + ylim(0,.8) +
    theme(legend.position = 'bottom',
          legend.title = element_text(size = 20), legend.text = element_text(size = 15),
          axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20),
          legend.spacing.x = unit(.7, 'cm'))+
    guides(color = guide_legend(override.aes = list(size = 2) ) )+
    scale_color_manual(values = c('Black','Red','Green','Blue'))+
    annotate('text',x=8.75,y=0.73,size=8,label = expression(paste(alpha,' = 1.95')))
