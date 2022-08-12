library(ggplot2)
library(gridExtra)
library(directlabels)
library(dplyr)
library(scales)
library(SPQR)
loc1 = load.SPQR('ffnn_45_1',path = 'EVP_local_1model')
loc3 = load.SPQR('ffnn_45_2',path = 'EVP_rmj')
loc1 = fit1
loc3 = fit2
loc2 = load.SPQR('ffnn_5',path = 'data/GP_global_ungridded')
p1 = plotGOF(loc1) + ggtitle(expression(paste('Q-Q plot - site 45 (',delta,' < 0.50)')))
p1
p2 = plotGOF(loc3) + ggtitle(expression(paste('Q-Q plot - site 45 (',delta,' > 0.50)')))
grid.arrange(p1,p2,nrow=1)
var.names <- c("delta","rho",paste0("X[",1:3,"]"))
p3 = plotQVI(loc3,var.index = 1:5,var.names = var.names)
# p3

b=p3$data
str(b)
tau = substr(b$tauexp,start = 6,stop=9)
tau = as.numeric(tau)
b$tau = tau
b = b %>% 
    mutate(x = recode_factor(x, `delta` = "delta", `rho` = "rho", `X[1]` = "X[1]", `X[2]` = "X[2]", `X[3]` = "X[3]"))
vi45_2 = ggplot(b,aes(x=tau,y=y,shape=x)) + geom_line(aes(col=x),size=1) +
    geom_point(size=3) + theme_bw() +
    labs(title = expression(paste('Variable importance - site 45 (',delta,' > 0.50)')),
         x = "Quantile",
         y = "Variable Importance",
         shape = "Covariate",
         col = "Covariate") + 
    scale_colour_discrete(labels = parse_format())+
    scale_shape_discrete(labels = parse_format()) + theme(legend.position = "bottom")

vi45_2
legend = cowplot::get_legend(vi45_2)
vi45_2 = ggplot(b,aes(x=tau,y=y,shape=x)) + geom_line(aes(col=x),size=1) +
    geom_point(size=3) + theme_bw() +
    labs(title = expression(paste('Variable importance at site 45 (',delta,' > 0.50)')),
         x = "Quantile",
         y = "Variable Importance") +
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))

vi45_2

p4 = plotQVI(loc1,var.index = 1:5,var.names = var.names)
# p3

b1=p4$data

str(b1)
tau = substr(b1$tauexp,start = 6,stop=9)
tau = as.numeric(tau)
b1$tau = tau
b1 = b1 %>% 
    mutate(x = recode_factor(x, `delta` = "delta", `rho` = "rho", `X[1]` = "X[1]", `X[2]` = "X[2]", `X[3]` = "X[3]"))
labels = c(expression(delta),expression(rho),expression(X[1]),expression(X[2]),expression(X[3]))

vi45 = ggplot(b1,aes(x=tau,y=y,shape=x)) + geom_line(aes(col=x),size=1) +
    geom_point(size=3) + theme_bw() + 
    labs(title = 'Variable importance at site 45',
         x = "Quantile",
         y = "Variable Importance",
         shape = "Covariate", col = "Covariate") +
    theme(legend.position = "right",plot.title = element_text(hjust = 0.5),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 14),
          title = element_text(size = 15))  +
    scale_color_manual(values = c('red','green','blue','magenta','black'),labels=labels)+
    scale_shape_manual(values = c(19,17,15,3,7),labels=labels)

vi45
grid.arrange(vi45,vi45_2,legend,nrow=2,ncol=2,
             layout_matrix = rbind(c(1,2),c(3,3)),
             widths = c(2.7,2.7), heights = c(2.6,.4))
