
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library("ggsci")
library(RColorBrewer)
library(cowplot)
library(ggalt)

######################################################################################## max derived allele count ######################################################################################

list.VariantType_1<-c(0,-0.0001,-0.001,-0.01)

data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s0=read.table(file="../data/fig3_data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s0")
data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s01=read.table(file="../data/fig3_data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s01")
data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s001=read.table(file="../data/fig3_data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s001")
data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s0001=read.table(file="../data/fig3_data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s0001")

fig3B_model1_r08_N1_s0_mdac<-ggplot(data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s0,aes(x=as.numeric(mdac)/2000,y=as.numeric(Mean),colour=factor(LD)))+geom_xspline(spline_shape=-0.6,aes(colour=factor(LD)))+geom_ribbon(aes(fill=factor(LD),ymin=(as.numeric(Mean)-as.numeric(Se)/1),ymax=(as.numeric(Mean)+as.numeric(Se)/1)),size=0,alpha=0.4)+theme_classic()+labs(x="",y="LD statitics (normalized, s=0)",title="B")+theme(legend.position = "bottom",legend.key.size=unit(10,"pt"),text=element_text(family="serif",size=12),axis.text=element_text(size=10,color="black"),legend.text=element_text(size=10),legend.title=element_text(size=10))+geom_point(size=0.6)+scale_y_continuous(expand=c(0,0),limits=c(0,8),breaks=seq(0,8,1))+geom_hline(yintercept=1,lty=2,size=0.8,color="gray")+scale_color_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))+scale_fill_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))
fig3B_model1_r08_N1_s01_mdac<-ggplot(data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s01,aes(x=as.numeric(mdac)/2000,y=as.numeric(Mean),colour=factor(LD)))+geom_xspline(spline_shape=-0.6,aes(colour=factor(LD)))+geom_ribbon(aes(fill=factor(LD),ymin=(as.numeric(Mean)-as.numeric(Se)/1),ymax=(as.numeric(Mean)+as.numeric(Se)/1)),size=0,alpha=0.4)+theme_classic()+labs(x="",y="LD statitics (normalized, s=0.01)",title="B")+theme(legend.position = "bottom",legend.key.size=unit(10,"pt"),text=element_text(family="serif",size=12),axis.text=element_text(size=10,color="black"),legend.text=element_text(size=10),legend.title=element_text(size=10))+geom_point(size=0.6)+scale_y_continuous(expand=c(0,0),limits=c(-1,4),breaks=seq(-1,4,0.5))+geom_hline(yintercept=1,lty=2,size=0.8,color="gray")+scale_color_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))+scale_fill_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))
fig3B_model1_r08_N1_s001_mdac<-ggplot(data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s001,aes(x=as.numeric(mdac)/2000,y=as.numeric(Mean),colour=factor(LD)))+geom_xspline(spline_shape=-0.6,aes(colour=factor(LD)))+geom_ribbon(aes(fill=factor(LD),ymin=(as.numeric(Mean)-as.numeric(Se)/1),ymax=(as.numeric(Mean)+as.numeric(Se)/1)),size=0,alpha=0.4)+theme_classic()+labs(x="",y="LD statitics (normalized, s=0.001)",title="B")+theme(legend.position = "bottom",legend.key.size=unit(10,"pt"),text=element_text(family="serif",size=12),axis.text=element_text(size=10,color="black"),legend.text=element_text(size=10),legend.title=element_text(size=10))+geom_point(size=0.6)+scale_y_continuous(expand=c(0,0),limits=c(-1,4),breaks=seq(-1,4,0.5))+geom_hline(yintercept=1,lty=2,size=0.8,color="gray")+scale_color_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))+scale_fill_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))
fig3B_model1_r08_N1_s0001_mdac<-ggplot(data.list.model1_h05_r08_diffMDac_N10000.samples1000.LD_s0001,aes(x=as.numeric(mdac)/2000,y=as.numeric(Mean),colour=factor(LD)))+geom_xspline(spline_shape=-0.6,aes(colour=factor(LD)))+geom_ribbon(aes(fill=factor(LD),ymin=(as.numeric(Mean)-as.numeric(Se)/1),ymax=(as.numeric(Mean)+as.numeric(Se)/1)),size=0,alpha=0.4)+theme_classic()+labs(x="",y="LD statitics (normalized, s=0.0001)",title="B")+theme(legend.position = "bottom",legend.key.size=unit(10,"pt"),text=element_text(family="serif",size=12),axis.text=element_text(size=10,color="black"),legend.text=element_text(size=10),legend.title=element_text(size=10))+geom_point(size=0.6)+scale_y_continuous(expand=c(0,0),limits=c(0,8),breaks=seq(0,8,1))+geom_hline(yintercept=1,lty=2,size=0.8,color="gray")+scale_color_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))+scale_fill_manual(values=c("#F8766D","#FA9FB5","#9ECAE1"))
plot_grid(fig3B_model1_r08_N1_s0_mdac,fig3B_model1_r08_N1_s0001_mdac,fig3B_model1_r08_N1_s001_mdac,fig3B_model1_r08_N1_s01_mdac,ncol=2,align="vh")

fig3B<-plot_grid(fig3B_model1_r08_N1_s0_mdac,fig3B_model1_r08_N1_s0001_mdac,fig3B_model1_r08_N1_s001_mdac,fig3B_model1_r08_N1_s01_mdac,ncol=2,align="vh")
ggsave("fig3B.maxdac.se.svg",fig3,width = 9, height = 6)

