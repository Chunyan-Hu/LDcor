#!/usr/bin/Rscript

## Rscript /sharedata/huchunyan/生物多样性信息学研究组Private/SLIs/scripts_calculation/001_5.test_V_vs_Va.single_func.r /sharedata/huchunyan/生物多样性信息学研究组Private/SLIs/basic_info/model1_h05_r05_N10000_1/ functype.S 5 model1_h05_r05 N10000.1run 2
## Rscript /media/sheldon/huchy/ani_from_zhou/20220816/001_5.test_V_vs_Va.single_func.r /media/sheldon/huchy/ani_from_zhou/2021.08.30/02.V_vs_Va.unrelated/panda/ vep 5 panda 2 

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
mindac <- args[5]
distance <- as.numeric(args[6]) # in bp

setwd(mpath)	
#names.vec <- c(0,-0.1,-0.01,-0.001,"-0.0001")
names.vec <- c("codingsynononly","missense","splicenonsense")

results <- matrix(ncol=11, nrow=length(names.vec))
results[,1]  = names.vec

for (funcii in c(1:length(names.vec))){
	## store load  and additive variance values for all dacs
	fname <- paste(mpath,group,".ALL.dac.frq.count.",names.vec[funcii],".",level, ".uppercaseaa.hist", corr, sep="")
	fname_position <-paste(mpath,group,".frq.count.indiv.",names.vec[funcii],".",level, ".uppercaseaa", corr, sep="")
	#print(paste("reading ",fname," ... ",proc.time(),sep=""))
	loads <- as.matrix(fread(fname,nThread=10))
	loads_position <- read.table(fname_position,fill=T,col.names=paste("V",1:(ncol(loads)+10)),sep="")[1:2];loads_position=as.matrix(loads_position)
	#print(loads_position[1:6,])
	colnames(loads) <- NULL # now we get a matrix with Nindivs rows, and each column is of a locus
	gnsample <- nrow(loads)
	loadsff <- loads[,colSums(abs(loads)) >= mindac] # we consider mutated loci only 
	loads_position <- loads_position[colSums(abs(loads)) >= mindac,]
	#print(loads_position[1:6,])
	rm(loads)
	loads <- loadsff[,colSums(abs(loadsff)) <= maxdac] # further, we can limit the maxdac
	loads_position <- loads_position[colSums(abs(loadsff)) <= maxdac,]
	rm(loadsff)
	#print(loads_position[1:6,])
	ref_chr=loads_position[1,1] ## 去除离得太近的snp
	ref_pos=loads_position[1,2]
	loci_to_remove=c(0)
	for (line in c(1:nrow(loads_position))){
		if (loads_position[line,1] == ref_chr){
			if (as.numeric(loads_position[line,2])-as.numeric(ref_pos) <= distance){loci_to_remove=c(loci_to_remove,line)}
			else {ref_chr=loads_position[line,1];ref_pos=loads_position[line,2]}
		}
	}
	loci_to_remove=loci_to_remove[-1]
	loads=loads[,-loci_to_remove]
		
	#print(paste(names.vec[funcii],": ",proc.time()," caculating...",sep=" -- "))
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
	#print(paste(names.vec[funcii],": ",proc.time()," caculation done",sep=" -- "))
	rm(loads)
}
outfile = paste(mpath,group,".maxdac",maxdac,".mindac",mindac,".samples", gnsample,".distance.",distance, sep="") # error
results.df <- data.frame(results)
colnames(results.df) <- c("VariantType", "Mean", "Va","V", "V/Va", "Nloci", "NetLD","cor","cor_abs","pos","neg")
#write.table(results.df, file = "", sep = "\t", row.names=FALSE, quote=FALSE)
write.table(results.df, file = outfile, sep = "\t", row.names=FALSE, quote=FALSE, append=FALSE)

