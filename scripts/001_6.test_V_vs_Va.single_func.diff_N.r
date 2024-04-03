#!/usr/bin/Rscript

## Rscript /media/sheldon/huchy/ani_from_zhou/20220608/scripts_calculation/001_6.test_V_vs_Va.single_func.diff_N.r /media/sheldon/huchy/ani_from_zhou/20220608/basic_info/model1_h00_r05_1/ functype.S 1 model1_h00_r05_1run 1 1000

library(boot)
library("matrixStats")
library("BSDA")
library(data.table)

####### read parameters ######
args <- commandArgs(trailingOnly = TRUE)
mpath <- args[1]
#funci <- args[2]
# splicenonsense
level <- args[2]
maxdac <- as.numeric(args[3])
#nsample <- as.numeric(args[3])
# 35
mpathnum=""
corr <- ""
group<- args[4]
##### set paramaters ############
#mindac <- 1
#names.vec<- args[5]
mindac <- args[5]
gnsample <- args[6] ## try 1000/100/500/50

setwd(mpath)
#names.vec <- c("-0.1", "-0.01", "-0.001","-0.0001","-1e-5","0")	
names.vec <- c(0,-0.1,-0.01,-0.001,"-0.0001","DFE")
results <- matrix(ncol=11, nrow=length(names.vec))
results[,1]  = names.vec

for (funcii in c(1:length(names.vec))){
	## store load  and additive variance values for all dacs
	fname <- paste(mpath,group,".",names.vec[funcii],".ALL.dac.frq.count.",level, ".uppercaseaa.hist", corr, sep="")
	print(paste("reading ",fname," ... ",proc.time(),sep=""))
	loads <- as.matrix(fread(fname,nThread=10))
	colnames(loads) <- NULL # now we get a matrix with Nindivs rows, and each column is of a locus
	##gnsample <- nrow(loads)
	loads <- loads[1:gnsample,] # we consider the header gnsample indivs only (try 1000/100/500/50)
	loadsff <- loads[,colSums(abs(loads)) >= mindac] # we consider mutated loci only and limit the mindac
	rm(loads)
	loads <- loadsff[,colSums(abs(loadsff)) <= maxdac] # further, we can limit the maxdac
	rm(loadsff)
	print(paste(names.vec[funcii],": ",proc.time()," caculating...",sep=" -- "))
	### compute statistics #####
	#print(loads)
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
	pos <- sum(Cor>0)-ncol(Cor)
	neg <- sum(Cor<0)	
	cor_abs <-(sum(abs(Cor))-ncol(Cor))/nloci/(nloci-1)

	results[funcii, 2] <- M
	results[funcii, 3] <- Va
	results[funcii, 4] <- V
	results[funcii,5] <- V/Va
	results[funcii,6] <- nloci
	results[funcii,7] <- Iscaled
	results[funcii,8] <- cor
	results[funcii,9] <- cor_abs
	results[funcii,10] <- pos
	results[funcii,11] <- neg
	print(paste(names.vec[funcii],": ",proc.time()," caculation done",sep=" -- "))
	rm(loads)
}
outfile = paste(mpath,group,".maxdac",maxdac,".mindac",mindac,".samples", gnsample, sep="") # error
results.df <- data.frame(results)
colnames(results.df) <- c("VariantType", "Mean", "Va","V", "V/Va", "Nloci", "NetLD","cor","cor_abs","pos","neg")
#write.table(results.df, file = "", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(results.df, file = outfile, sep = "\t", row.names=FALSE, quote=FALSE, append=FALSE)
	
