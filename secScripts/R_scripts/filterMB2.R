#!/usr/bin/env Rscript

#script to post filter genes in MGS to identify "core" genes
#(c) Falk Hildebrand

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

if (FALSE){
	args=c()
	args[1] = "C:/Users/hildebra/OneDrive/science/data/test/Cano_Alien//SB.clusters"
}


inP = args[1]
inObs = paste0(inP,".obs")
inObs=gsub("Wclusters","clusters",inObs)
print(inObs)

picD=paste0(dirname(inP),"/pics/")
dir.create(picD, showWarnings = FALSE)
if (length(args)<2) {#output
	args[2] = paste0(args[1],".core")
} else {outP = args[2]}
if (length(args)<3) {
	args[3] = paste0(args[1],".ext")
} 
genesUsed = list()

print(args[1])
M=read.table(args[1],TRUE,"\t",as.is=TRUE)
#all important criteria for which  gene is acceptable at all..
coreCriteria = (M[,4]/M[,3] < 0.1 & M[,5] > 0.9 & M[,5] < 3 ) | M[,6]

print(M[1,])

print(paste0(sum(coreCriteria),"/",dim(M)[1]," Genes passed Core Criteria"))



Meta = as.matrix(read.table(inObs,FALSE,"\t",as.is=TRUE))
Obs = as.numeric(Meta[,2]);names(Obs)=Meta[,1]
Obs[is.na(Obs)]=0
Mret=matrix(NA,3,0)
bins=unique(M[,1]);b=bins[2] # 37 45 48
minObs=Obs


pdf(paste0(picD,"BinPlots.pdf"),7,7)#740,740)

for (b in bins){
	sel=M[,1]%in%b & coreCriteria
	maxv = Obs[b]
	vals=M[sel,3]
	isMG = M[sel,6]
	hiMG = quantile(vals[isMG],prob=c(0.95))
	#if (length(vals)>20000){vals=vals[1:10000]}
	mv=quantile(vals,prob=c(0.97))
	#maxv = vals[10]
	#if (maxv>vals[10]){maxv=vals[10]}
	if (maxv>mv){maxv=mv}
	minObs[b] = max(maxv*0.1,2)
	corecut = quantile(vals,prob=c(0.5,0.99))
	stidx = sum(vals>maxv)
	if (maxv<10){
		retidx =vals >= corecut[1] & vals <= corecut[2] #%in% unique(vals[stidx:min((1000+stidx),length(vals)*0.7)])
		if (maxv<4){minObs[b]=1}
		#plot((vals));abline(h=corecut)
		#	IQR(vals)
	} else {
		#blast= b #DEBUG
		hv=hist(vals,breaks =40,plot=FALSE)
		#hv=hist(vals,breaks =40,plot=TRUE)
		#str(hv)
		hiHists = hv$breaks[which(hv$counts>quantile(hv$counts,probs=c(0.6)))]
		hiSel = hiHists >= maxv*0.5
		lowVal = max(2,maxv*0.1)
		if (sum(hiSel)<2){
			#vals >= corecut[1] & vals <= corecut[2]#
			retidx = vals %in% unique(vals[stidx:min((1000+stidx),length(vals)*0.7)])
		} else {		
			rhi=round(range(hiHists[hiSel]))
			#coreVals = vals >= corecut[1] & vals <= corecut[2]#
			#length(coreVals)
			#coreVals = vals[ vals>rhi[1] & vals<rhi[2]] #
			#rhi2 = c(median(coreVals)-2*IQR(coreVals),median(coreVals)+2*IQR(coreVals))
			rhi2=rhi
			rhi2[2]= max( rhi2[2] + ((maxv-rhi2[2])/2)       ,  hiMG)
			#rhi2[2] = maxv
			#if (rhi2[2] > rhi[2]) {rhi2[2] = rhi[2]}
			retidx=vals>=rhi2[1] & vals<=rhi2[2]
			#retidx = vals >= corecut[1] & vals <= corecut[2]
			#plot(density(coreVals))
			plot((vals),main=paste(b," (",paste(rhi2,collapse=";"),")"), col = isMG+1,log="x",pch=as.numeric(retidx)+2);
			grid(col="gray",lty=2)
			#abline(h=rhi,col="darkgreen");
			abline(h=rhi2,col="darkred");#abline(h=lowVal,col="darkblue")
			#abline(h=corecut[1],col="green");abline(h=corecut[2],col="red");
			abline(h=maxv,col="darkred",lty=2,lwd=2)
			text(sum(retidx),rhi2[1],sum(retidx),adj=c(-1,1))
			text(length(vals)*0.98,rhi2[1],length(vals),adj=c(1,0))
		}
		
		#library(mixtools)
		#nem=spEMsymloc(vals,mu=c(2,370))
		#numNs = length(nem$mu)
		#denV=density(vals[vals>mv],adj=1.7)
		#nem=mvnormalmixEM(rbind(denV$x,denV$y),k=2)
		#max=denV$x[which.max(denV$y)]
		#plot(denV)
		#abline(v=nem$mu,col="green")
		#abline(v=nem$mu,col="green")
		#points(nem$mu+nem$sigma,c(0,0))
	}
	retidx = retidx# | isMG
	#M[which(sel)[retidx],]
	genesUsed[[b]] = sum(retidx)
	Mret=rbind(Mret, M[which(sel)[retidx],] )
	#hist(vals,plot=FALSE,breaks=10)
	#quantile(vals,prob=c(0.5,0.9))
}

dev.off()

pdf(paste0(picD,"NumGenesPerBinCore.pdf"),12,7)#1500,740)
plot(unlist(genesUsed),type="p")
dev.off()

#extended core genome - basically everything I would consider truly part of a species
extCriteria = (M[,5] > 0.9 & M[,3] >= minObs[M[,1]]) | M[,6]
Mext = M[extCriteria,]
pdf(paste0(picD,"NumGenesPerBinExt.pdf"),12,7)#1500,740)
plot(table(Mext[,1]),type="p")
dev.off()


# last: tables
write.table(Mret,args[2],quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
write.table(unlist(genesUsed),paste0(args[2],".cnts"),quote=FALSE,sep="\t",row.names=TRUE,col.names=FALSE)
write.table(Mext,args[3],quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
