#script takes mitag tables from MATAFILER and normalizes by metagStat.txt
#(c) Falk Hildebrand
args = commandArgs(trailingOnly=TRUE)
inD = args[1]
inD=paste0(inD,"/")
meta=read.table(paste0(inD,"metagStats.txt"),sep="\t",quote="",header=1,row.names=1,comment.char="",as.is=1)

rdSum = meta[,"Accepted1"]+meta[,"Accepted2"]
names(rdSum) = dimnames(meta)[[1]]

SUtag=c("SSU","LSU")
lvl=c("domain","phylum","class","order","family","genus","species","hit2db")
for (y1 in 1:length(SUtag)){
	for(yy in 1:length(lvl)){
		inMF = paste0(inD,"/pseudoGC/Phylo/RiboFind/",SUtag[y1],".miTag.",lvl[yy],".txt.gz")
		savF = paste0(inD,"/pseudoGC/Phylo/RiboFind/",SUtag[y1],".miTag.",lvl[yy],".Rdata")
		if (!file.exists(inMF)){
			print (paste("Can't find ",inMF))
			next
		}
		M=read.table(
			gzfile(inMF)
			,header=TRUE,sep="\t",quote="",comment.char="",row.names=1,as.is=TRUE)
		dimnames(M)[[2]] = gsub("\\.[^\\.]+\\.hiera\\.txt\\.gz","",dimnames(M)[[2]])
		rdSum2=rdSum[dimnames(M) [[2]] ]; rdSum2[is.na(rdSum2)] = 0
		subs=colSums(M)>0 & rdSum2>0
		Mt = t(t(M[,subs])/rdSum2[subs])
		save(Mt,rdSum,file=savF)
	}
}