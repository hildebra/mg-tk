if(!require("ape",quietly=TRUE,warn.conflicts =FALSE)){install.packages("ape",repos="https://cloud.r-project.org");require("ape")}
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#inF= "C:/Users/falkh/OneDrive/science/data/test/tree/RXng_allsit.raxml.bestTree";tar="specI_v2_0150"
inF = args[1]
tar = args[2]
tree=read.tree(inF)
cop = cophenetic(tree)
c1 = cop[tar,]
sc1=sort(c1)
res=sc1[sc1>0.01][1:10]
cat(paste0(paste(names(res),collapse=" "),"\n"))





