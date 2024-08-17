#!/usr/bin/env Rscript

#all packages required to run MATAFILER & extended options, execute to install (and don't conflict later)

#most of these are already installed via conda/mamba

#for remainder, sometimes there's a problem with the path


args = commandArgs(trailingOnly=TRUE)
customRlibP = NULL
if (length(args)>0){
	customRlibP = args[1]
	print(paste("Using custom R lib path:",customRlibP))
}


libDirs = .libPaths()
instLib = libDirs[grep("MGTK",libDirs)] #set preferential dir to MFF env..

print (instLib)

if (file.access(instLib, mode = 2) == -1){
	for (xy in .libPaths()){
		instLib = xy
		print (xy)
		if (file.access(instLib, mode = 2) == 0){
			print (paste("Can't access file path",instLib))
			break
		}
	}
}

if (file.access(instLib, mode = 2) == -1){
	stop(paste("Could not find writeable path in " , paste(.libPaths(), collapse=", ")))
}

#if (length(grep("MFF",.libPaths())) >= 1){
	#.libPaths( c(.libPaths()[2], .libPaths()[1]) )
	#.libPaths(.libPaths()[grep("MFF",.libPaths())] )
#	instLib = .libPaths()[grep("MFF",.libPaths())] 
#}
print(paste0("Installing R libraries to dir: ",instLib))

suppressPackageStartupMessages({

if(!require("argparser",quietly=TRUE,warn.conflicts =FALSE)){install.packages("argparser",repos="https://cloud.r-project.org");require("compiler")}
if(!require("optparse",quietly=TRUE,warn.conflicts =FALSE)){install.packages("optparse",lib=instLib,repos="https://cloud.r-project.org");require("optparse",warn.conflicts=FALSE, quietly = TRUE)}
#if(!require("compiler",quietly=TRUE,warn.conflicts =FALSE)){install.packages("compiler",lib=instLib, repos="https://cloud.r-project.org");require("compiler")}


if(!require("dplyr",quietly=TRUE,warn.conflicts =FALSE)){install.packages("dplyr",lib=instLib, repos="https://cloud.r-project.org");require("dplyr",warn.conflicts=FALSE, quietly = TRUE)}
if(!require("RColorBrewer",quietly=TRUE,warn.conflicts =FALSE)){install.packages("RColorBrewer",lib=instLib, repos="https://cloud.r-project.org");require("RColorBrewer",warn.conflicts=FALSE, quietly = TRUE)}


if(!require("ape",quietly=TRUE,warn.conflicts =FALSE)){install.packages("ape",lib=instLib, repos="https://cloud.r-project.org");require("ape",warn.conflicts=FALSE, quietly = TRUE)}




if (!require("BiocManager", warn.conflicts=FALSE,quietly = TRUE)){install.packages("BiocManager",lib=instLib, repos="https://cloud.r-project.org");require("BiocManager");}
if (!require("ggtree", warn.conflicts=FALSE,quietly = TRUE)){require("BiocManager");BiocManager::install("ggtree");}

if(!require("phytools",warn.conflicts=FALSE,quietly=TRUE)){install.packages("phytools",lib=instLib, repos="https://cloud.r-project.org");require("phytools")}

if(!require("tidytree",warn.conflicts=FALSE,quietly=TRUE)){install.packages("tidytree",lib=instLib, repos="https://cloud.r-project.org");require("tidytree")}

if(!require("phangorn",quietly=TRUE,warn.conflicts =FALSE)){install.packages("phangorn",lib=instLib, repos="https://cloud.r-project.org");require("phangorn",warn.conflicts=FALSE, quietly = TRUE)}


}) #suppressPackageStartupMessages