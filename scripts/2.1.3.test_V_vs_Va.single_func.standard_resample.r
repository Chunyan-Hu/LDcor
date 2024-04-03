#!/usr/bin/Rscript

## Rscript /media/sheldon/huchy/ani_from_zhou/2021.08.30/2.1.3.test_V_vs_Va.single_func.standard_resample.r /media/sheldon/huchy/ani_from_zhou/2021.08.30/02.V_vs_Va.unrelated/panda/ vep 2 panda 2 splicenonsense missense 100

library(boot)
library("matrixStats")
library("BSDA")
library(data.table)

####### read parameters ######
args <- commandArgs(trailingOnly = TRUE)
mpath <- args[1] # /media/sheldon/huchy/ani_from_zhou/2021.08.30/02_1.V_vs_Va_with_kinship/
#funci <- args[2]
# splicenonsense
level <- args[2] # vep
maxdac <- as.numeric(args[3]) # 2
#nsample <- as.numeric(args[3])
# 35
mpathnum=""
corr <- ""
group<- args[4] # panda.v3/
##### set paramaters ############
mindac <- args[5] # 2
names.vec<- args[6]  # splicenonsense
pool.vec<- args[7] # codingsynononly
nround <- as.numeric(args[8]) # 100
#the_tail <- args[9]
setwd(mpath)
results <- matrix(ncol=14, nrow=nround+1)
results[,1]  = pool.vec
results[1,1]  = names.vec

fname <- paste(mpath,group,".ALL.dac.frq.count.", names.vec,".",level, ".uppercaseaa.hist",corr, sep="")
loads <- as.matrix(fread(fname))
#print(fname)
colnames(loads) <- NULL
loadsff <- loads[,colSums(abs(loads)) >= mindac] ## attention
rm(loads)
loads <- loadsff[,colSums(abs(loadsff)) <= maxdac]
rm(loadsff)
datavec <- rowSums(loads)
temp <- colSums(loads)
gnsample <- nrow(loads)
gnloci <- ncol(loads)
datavec_count <-temp
V <- var(datavec)
M <- mean(datavec)
Va <- sum(colVars(loads))
I <- V - Va
nloci <- sum(temp>0)
if (nloci > 1) Iscaled <- I/nloci/(nloci-1) else Iscaled<-NA
Cor <- cor(loads)
cor <-(sum(Cor)-ncol(Cor))/nloci/(nloci-1)
pos <- sum(Cor>0)-ncol(Cor)
neg <- sum(Cor<0)
cor_abs <-(sum(abs(Cor))-ncol(Cor))/nloci/(nloci-1)
results[1, 2] <- M
results[1, 3] <- Va
results[1, 4] <- V
results[1,5] <- V/Va
results[1,6] <- nloci
results[1,7] <- Iscaled
results[1,8] <- cor
results[1,9] <- cor_abs
results[1,10] <- pos
results[1,11] <- neg
results[1,12] <- 0
results[1,13] <- 0
results[1,14] <- 0
thecor=cor
thecorabs=cor_abs
#print(results[1,])
frq_dis <- c(mindac:maxdac)
for (ii in mindac:maxdac){frq_dis[ii]<-length(which(temp==ii))}
#print(frq_dis)

## store load  and additive variance values for all dacs
fname <- paste(mpath,group,".ALL.dac.frq.count.", pool.vec,".",level, ".uppercaseaa.hist", corr, sep="")
#print(paste("reading ",fname," ... ",proc.time(),sep=""))
loads <- as.matrix(fread(fname,nThread=10))
#print(fname)
colnames(loads) <- NULL
loadsff <- loads[,colSums(abs(loads)) != 0]
rm(loads)
loads0 <- loadsff[,colSums(abs(loadsff)) <= maxdac]
rm(loadsff)
temp <- colSums(loads0)
count_ID <- vector("list",maxdac)
for (ii in mindac:maxdac){count_ID[[ii]]<-which(temp==ii);print(length(count_ID[[ii]]))}
#Nround <- matrix(ncol=maxdac, nrow=nround)
for (funcii in c(1:nround)){
	
	## resample by distribution of syn
	#print(paste("round",funcii,": ",proc.time()," caculating...",sep=" "))
	loads <- matrix(data=NA, nrow=gnsample)	
	for (ii in mindac:maxdac){
		## Nround[funcii,ii] <- sample(count_ID[[ii]],frq_dis[ii]) # error
		loads <- cbind(loads,loads0[,sample(count_ID[[ii]],frq_dis[ii])])
	}
	loads <- loads[,-1]

	### compute statistics #####
	##print(loads)
	datavec <- rowSums(loads)
	temp <- colSums(loads)
	  
	## calculate variance and mean
	V <- var(datavec)
	M <- mean(datavec)
	  
	##calculate additive variance
	Va <- sum(colVars(loads))
	#Vloga <- sum(log10(colVars(loads)))
	  
	## calculate Net LD
	I <- V - Va
	nloci <- sum(temp>0)
	if (nloci > 1) Iscaled <- I/nloci/(nloci-1)
	else Iscaled<-NA
	
	Cor <- cor(loads)
	cor <-(sum(Cor)-ncol(Cor))/nloci/(nloci-1)
	cor_abs <-(sum(abs(Cor))-ncol(Cor))/nloci/(nloci-1)
	pos <- sum(Cor>0)-ncol(Cor)
	neg <- sum(Cor<0)	

	results[funcii+1, 2] <- M
	results[funcii+1, 3] <- Va
	results[funcii+1, 4] <- V
	results[funcii+1,5] <- V/Va
	results[funcii+1,6] <- nloci
	results[funcii+1,7] <- Iscaled
	results[funcii+1,8] <- cor
	results[funcii+1,9] <- cor_abs
	results[funcii+1,10] <- pos
	results[funcii+1,11] <- neg
	results[funcii+1,12] <- (thecor-cor)/abs(cor) 
	results[funcii+1,13] <- (thecorabs-cor_abs)/cor_abs
	results[funcii+1,14] <- (thecor-cor)/cor_abs
	#print(paste(funcii,": ",proc.time()," caculation done",sep=" -- "))
	#print(results[funcii+1,])
	rm(loads)
}
outfile = paste(mpath,group,".",names.vec,"_resampled_from_",pool.vec,".maxdac",maxdac,".mindac",mindac,".nrounds", nround,sep="")
results.df <- data.frame(results)
colnames(results.df) <- c("VariantType", "Mean", "Va","V", "V/Va", "Nloci", "NetLD","cor","cor_abs","pos","neg","sli1","sli2","sli3")
#write.table(results.df[1:11,], file = "", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(results.df, file = outfile, sep = "\t", row.names=FALSE, quote=FALSE)
	
