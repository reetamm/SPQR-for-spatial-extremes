library(gridExtra)
library(ggplot2)
loc   <- create_locations()
loc11 = loc$nn_id[11,]
loc18 = loc$nn_id[18,]
loc45 = loc$nn_id[45,]
loc18
loc80
n11 = rep(0,50)
n11[11] = 1
n11[loc11] = 2

n18 = rep(0,50)
n18[18] = 1
n18[loc18] = 2

n45 = rep(0,50)
n45[45] = 1
n45[loc45] = 2

locations = data.frame(loc$locs)
names(locations)
locations$n18 = factor(n18)
locations$n11 = factor(n11)
locations$n45 = factor(n45)
locations$labels = 1:50
p1 = ggplot(locations,aes(x=X1,y=X2,col=n11,shape=n11,label=labels)) + geom_point(size=2)+  
    scale_shape_manual(values=c(1,15,16)) + 
    scale_color_manual(values=c("black","blue", "red")) +
    labs(x="",y="") +
    theme_bw() + theme(legend.position = "none")+ geom_text(hjust=0, vjust=0)
p2 = ggplot(locations,aes(x=X1,y=X2,col=n18, shape=n18, label = labels)) + geom_point(size=2)+  
    scale_shape_manual(values=c(1,15,16)) + 
    scale_color_manual(values=c("black","blue", "red")) +
    labs(x="",y="") +
    theme_bw() + theme(legend.position = "none") #+ geom_text(hjust=0, vjust=0)

p3 = ggplot(locations,aes(x=X1,y=X2,col=n45, shape=n45, label = labels)) + geom_point(size=2)+  
    scale_shape_manual(values=c(1,15,16)) + 
    scale_color_manual(values=c("black","blue", "red")) +
    labs(x="",y="") +
    theme_bw() + theme(legend.position = "none") #+ geom_text(hjust=0, vjust=0)
p2


grid.arrange(p2,p3,ncol=2)

