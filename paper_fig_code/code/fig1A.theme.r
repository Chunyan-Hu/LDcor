
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library("ggsci")
library(RColorBrewer) 
library("ggalt")


data1<-data.frame(rbind(cbind(number=c(0:100),rep("LDcor=0"),c(dpois(0:100,5)/sum(dpois(0:100,5)))),cbind(number=c(0:100),rep("LDcor<0"),c(dpois(0:9,5)/sum(dpois(0:9,5),dpois(10:100,5)*0.95),dpois(10:100,5)*0.95)/sum(dpois(0:9,5),dpois(10:100,5)*0.95)),cbind(number=c(0:100),rep("LDcor>0"),c(dpois(0:9,5)/sum(dpois(0:9,5),dpois(10:100,5)*1.05),dpois(10:100,5)*1.05)/sum(dpois(0:9,5),dpois(10:100,5)*1.05))))

data1_1=data1[which(data1$V2=="LDcor=0"),]
data1_2=data1[which(data1$V2=="LDcor>0"),];data1_2$V3=as.numeric(data1_2$V3)-as.numeric(data1_1$V3)
data1_3=data1[which(data1$V2=="LDcor<0"),];data1_3$V3=as.numeric(data1_3$V3)-as.numeric(data1_1$V3)
data1_1$V3=as.numeric(data1_1$V3)-as.numeric(data1_1$V3)
data1_0=rbind(data1_1,data1_2,data1_3)
options(scipen=200)

figS14<-ggplot(data1_0,aes(x=as.numeric(number),y=as.numeric(V3),color=V2,fill=V2))+geom_point(size=0.8)+geom_xspline(size=0.7,spline_shape=-0.3)+theme_classic()+theme(legend.position=c(0.8,0.8),legend.key.size=unit(10,"pt"),text=element_text(family="serif",size=12),axis.text=element_text(size=10,color="black"),legend.text=element_text(size=10),legend.title=element_text(size=10))+labs(x="Mutation number",y="Frequency")+scale_y_continuous(expand=c(0,0),limits=c(-0.001,0.001),breaks=seq(-0.001,0.001,0.0005))+scale_x_continuous(limits=c(0,18),breaks=seq(0,18,1))+scale_color_manual(values=c(brewer.pal(6, "Blues")[3],brewer.pal(5,"Greys")[3],brewer.pal(9,"RdPu")[4]))+scale_fill_manual(values=c(brewer.pal(6, "Blues")[3],brewer.pal(5,"Greys")[3],brewer.pal(9,"RdPu")[4]))+geom_bar(stat="summary",fun=mean,width=0.2,alpha=0.6,size=0)
ggsave("fig1A.theme.svg",figS14,width = 9, height = 4)


