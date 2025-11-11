




# MG-TK manual

MG-TK uses three primary phases to analyse metagenomes: 
	1. initial raw sequence processing via MG-TK.pl
	2. building a genecatalog via secScripts/geneCat.pl
	3. reconstructing species from MAGs and inferring intraspecific phylogenies via secScripts/MGS.pl
	4. (building phylogenies via secScripts/phylo/buildTree5.pl)

### MG-TK.pl flags

<details>
  <summary>Expand section</summary>

```{sh}
# base flags
	-map $MAP 					mapping file, that also contains input and output directories
	-config $CONFIG			alternative config file  (Default: autodetect)

# flow related
	-submit [0|1]				0: dry run to check file paths, general submission works; 1: submit actual jobs, this will take a long time in most cases
	-from [#]					run subsample of mapping file starting at sample # (use with -to)
	-to [#]						run subsample of mapping file ending at sample # (use with -from)
	-ignoreSmpls [string]			comma separated list of #SmplIDs that are skipped (sample id in .map file)
	-rmSmplLocks [0|1]				1: remove existing sample locks (useful if jobs have crashed, leaving abondened sample locks) 
	-redoFails [0|1]				if any step of requested analysis failed, just redo everything (use with care!) 
	-maxConcurrentJobs [#]			max jobs in queue, useful for large samples sets, currently only works on slurm (see also -killDepNever)
	-killDepNever [0|1]				kills jobs in the "JobDependencyNeverMet" state, as these will block [maxConcurrentJobs], 	
	-excludeNodes [string]					exclude certain HPC nodes, comma separated list e.g. node1,node2,..
	-submSystem [qsub,SGE,bsub,LSF]	set submission system (default: autodetect)
	-redoContigStats [0|1]				if any step of requested analysis failed, contigStats (coverage per gene, kmers, GC content) will be deleted & started again
	-loopTillComplete [X:Y]			script will loop over the assigned samples until all jobs are finished #use synatx "X:Y" where X is num loops, Y is the window size, eg "6:250" would run 6 loops of max 250 samples, then move on to next 250 samples (#dangerous flag)
	-requireInput [0/1]		in case input reads are not present (e.g. something wrong in map), 0 will continue pipeline, 1 will abort
	-silent [0/1] 			Controls how much information is printed on console
	-OKtoRWassGrps [0|1]			1: can delete assemblies, if suspects error in them, powerful, but careful! (Default: 0)
	-maxUnzpJobs [#]				#how many unzip jobs to run in parallel (not to overload HPC IO). Default:20
	-skipSmallSmplsMB [#]			skip sample if the overall file size is < than given number (in MB). Default: 1
	-forceWriteStats [0|1]			force (re)writing of the metagStats HTML report and text files. Default: 0


# Detecting raw input files
	-inputFQregex1 [‘R1’]			R1 input regex extension (e.g. R1 could be '.*_1\.f[^\.]*q\.gz$' or last resort '(.*_pe_1\.f[^\.]*q\.gz$)|(.*R1_00\d\.f[^\.]*q\.gz$)|(.*[\._]1\.f[^\.]*q\.gz$)|(.*R1\.fq\.gz)' )
	-inputFQregex2 [‘R2’]			R2 input regex extension (e.g. R2 could be '.*_2\.f[^\.]*q\.gz$' or last resort '(.*_pe_2\.f[^\.]*q\.gz$)|(.*R2_00\d\.f[^\.]*q\.gz$)|(.*[_\.]2\.f[^\.]*q\.gz$)|(.*R2\.fq\.gz)' )
	-inputFQregexSingle[‘RS’]			RS Single read regex extension, see -inputFQregex1 for examples
	-inputFQregexTrustSingle [0/1]		(1) if grep of files (rawSrchString) has multi assignments, the -inputFQregexSingle takes precendence over -inputFQregex1. (Default: 0)
	-inputBAMregex[‘BS’]				BS regex to detect Bam files (e.g. '.*\.bam$'). Currently only implemented for unpaired reads (eg PacBio output)

# file structure
	-rm_tmpdir_reads [0|1]			remove tmpdir with reads (default: 1)
	-rm_tmpInput [0|1]			remove raw, human / adaptor filtered reads, if sdm clean created? (and not needed any longer)
	-reduceScratchUse 1			remove scratch dir; should always be 1, unless debugging
	-globalTmpDir $PATH			absolute path to global shared tmp dir (like a scratch dir)
	-nodeTmpDir $PATH			absolute path to tmp dir on local HDD of each executing node
	-nodeHDDspace $PATH			HDD tmp space to be requested for each node. Some systems don't support this
	-legacyFolders 0				output folders use read dir as name (1) or mapping file (0). default: 0

# preprocessing (cleaning reads etc/input FQ related)
	-useTrimomatic [0|1]			remove adapter seq from input reads (default: 0 as sdm does this)
	-usePorechop [0|1]			adapter removal for Nanopore, should be automatic activated
	-splitFastaInput [0|1]			split fasta input
	-mergeReads [0|1] 			merge reads
	-ProbRdFilter [???]				sdm probabilistic filter ???
	-pairedReadInput [0|1]			0: not paired input 1: read pairs are expected in each in dir
	-inputReadLength [#]			read length #
	-filterHostRds				0: do not filter host reads, 1: filter host reads via kraken2, 2: .. via kraken1, 3: .. via hostile
	-filterHostKrak2DB $PATH		path to host kraken database
	-hostileIndex $PATH			name of hostile DB to filter against (default: human-t2t-hla)
	-onlyFilterZip [0|1]			??
	-mocatFiltered [0|1]			??
	-filterHostKr2Conf [#]                      set host kraken2 confidence parameter (e.g. 0.05) 
	-filterHostKr2Quick ["--quick"]             set quick option for kraken2 (should be "--quick")
	-upload2EBI [dir]			copy reads after host read removal to this directory, ready for upload to standard repositories

# sdm (read filtering) related
	-gzipSDMout [0|1] 			gzip sdm output
	-XfirstReads [int]				only use X first reads of each input read file
	-minReadLength [#]			minimum read length in sdm filtering step
	-maxReadLength [#]			maximum read length in sdm filtering step

# assembly
	-assembleMG [1|2|3|4|5] 			which assembler to use; 1: spades; 2: Megahit 3: FLYE 4: metaMDBG 5: hybrid assemblies (megahit+metaMDBG)
	-assemblCores [#] 				number of cores # for assembly (same as -spadesCores); e.g. 12
	-assemblyKmers [“#1,#2,#3,#4,…”]	number of kmers for assembly, comma-delimited e.g. "25,43,67,87,111,131" 
	-assemblMemory [#]			memory used for assembly (same as -spadesMemory); e.g. 100
	-asssemblyHddSpace [#]			HDD space requested by assembler in Gb; e.g. 120. Default: auto

# gene prediction on assembly
	-predictEukGenes [0|1]		predict eukaryotic genes; severely limits the total predicted gene amount (~25% of total genes) (Default: 0)
	-kmerPerGene [0|1]			1: report kmer frequencies per gene (Default: 0)

# mapping
	-mapper [1|2|3|4]				1: bowtie2, 2:bwa, 3: minimap2, 4:kma, 5:strobealign -1:auto (bowtie2 short, minimap2 long reads), -2:auto(strobealign short, minimap2 long). (Default: -1)
	-mappingCores [#]				cores # used for mapping
	-mappingMem [#]				memory # used for mapping bwa/bwt2 in GB (Default: auto (-1))
	-mapSortMem  [#]			memory # used for samtools sort in GB (Default: auto (-1))
	-mappingCoverage			1: calculate coverage per predicted gene, contig, windows (Default: 1)
	-mapSupportReadsOntoAssembly [#]   1: map also support reads (e.g. PacBio in case of hybrid assemblies) onto final assembly. (Default: 0)
	-rmDuplicates [0|1]			1: remove read duplicates (Default: 1)
	-mapperFilterIll [# # #]				parameters for postprocessing mappings for short read data (Default: "0.05 0.75 20")
	-mapReadsOntoAssembly [0|1] 	1: map original reads to assembly to estimate contig (and bin, gene etc) abundance (Default: 1)
	-saveReadsNotMap2Assembly [0|1]			1: save reads not mapping to assembly in separate file (Default: 0)
	-remap2assembly [0|1]			1: redo the mapping to assembly (Default: 0)
	-JGIdepths [0|1]			1: calculate jgi coverage, only required when using MetaBAT2 binning (Default: automatic)
	-mapperLargeRef [0|1]			1: reference DB that is mapped against is veryyy large (highly unusual that this is needed) (Default: 0)
	-mapSaveCRAM [0|1]				1: keep .cram from mapping reads to assembly, 0: delete after essential processes have used them (Default: 1)

# binning
	-Binner [1|2|3]				#0=no binning, 1= do metaBat2 binning, 2=SemiBin, 3= MetaDecoder (experimental) (Default: 0)
	-BinnerMem [#]			define binning memory; e.g. 600 (Default: automatic)
	-BinnerCores [#]			cores used for Binning process (and checkM)
	-redoEmptyBins [0|1]		mostly for debugging: redo every sample where no bins where found (note that in some metagenomes there might be no bins)
	-checkM2 [0|1]			using checkM2 to assess bin quality
	-checkM1 [0|1]			using checkM1 to assess bin quality
	-redoBinning [0|1]			redo binning

# SNPs
	-getAssemblConsSNP [0|1]		1: SNPs (onto self assembly); calculates consensus SNP of assembly (useful for checking assembly gets consensus and Assmbl_grps)
	-get2ndMappingConsSNP [0|1]		1: calculate consensus SNPs for mappings against references (map2tar mode)
	-redoAssmblConsSNP [0|1]		1: redo getAssemblConsSNP (Default: 0)
	-redoGeneExtrSNP [0\1]			1: redo gene extractions from consensus SNP contig (Default: 0)
	-SNPjobSsplit [#]				#: split consensus SNP jib further (Default: 1)
	-SNPsaveVCF [0|1]				1: save SNPs to VCF file (Default: 1)
	-SNPsaveConsFasta [0|1]			1: save consensus fasta based on SNP calls (Default: 0)
	-SNPcaller [MPI|FB]				"MPI" mpileup or ".FB" for freebayes (Default: MPI)
	-SNPcores [#]				number of cores ‘#’ used for SNP calling
	-SNPmem [#]				Memory allocated for consensus SNP calling process in Gb (Default: 23)
	-SNPconsMinDepth [#]			how many reads coverage to include position for consensus call? (Default: 0)
	-SNPminCallQual [#]				min Quality of a SNP call to be accepted (Default: 20)
	-SVcaller [0|1]					call structural variants in assembly. 1=delly, 2=gridss. (Default: 0).

# functional profiling (raw reads without assemblies)
	-profileFunct [0|1]				1: do diamond functional profiling (Default: 0)
	-reParseFunct [0|1]			1: redo diamond result parsing and translation to categories (Default: 0)
	-reProfileFunct [0|1]			redo diamond functional profiling  (Default: 0)
	-reProfileFuncTogether [0|1]		if any func database needs to be redone, than redo all indicated databases (useful if number of reads changed)  (Default: 0)
	-diamondCores [#]				number of cores ‘#’ used for diamond run (Default: 12)
	-DiaParseEvals [#]				evalues at which to accept hits to func database (Default: 1e-7)
	-DiaSensitiveMode [0|1]			1: run diamond in sensitive mode  (Default: 0)
	-rmRawDiamondHits [0|1]		1: remove raw diamond hits  (Default: 0)
	-DiaMinAlignLen [#]			set diamond min alignment length to accept hit ‘#’  (Default: 20)
	-DiaMinFracQueryCov [f]			set diamond minimum fraction query coverage ‘f’  (Default: 0.1)
	-DiaPercID [#]				set diamond percent identify ‘#’ to accept hit  (Default: 40)
	-diamondDBs [db]				set diamond ‘db’ using one or more of: NOG,MOH,ABR,ABRc,ACL,KGM,CZy,PTV,PAB,MOH2 , can be comma separated for multiple DBs

# ribo profiling (miTag)
	-profileRibosome [0|1]			1: profile ribosomal SUs, (Default: 0)
	-riobsomalAssembly [0|1]		1: assemble ribosomal SUs, (Default: 0)
	-reProfileRibosome [0|1]			1: redo -profileRibosome, (Default: 0)
	-reRibosomeLCA [0|1]			1: redo classify ribosomal SUs (Default: 0)
	-riboMaxRds [#]				Number of presorted reads to use in LCA to control computational time (Default: 50,000)
	-saveRiboRds [0|1]			1: store raw presorted hits (Default: 0)
	-thoroughCheckRiboFinish [0|1]	1: check ribo profiling was successful (Default: 0)

# other tax profilers..
	-profileMetaphlan2 [0|1|3]		3: perform metaphlan3 (check local metaPhlan version) read profiling (Default: 0)
	-profileMOTU2 [0|1]			1: perform mOTUs2/mOTUs3 (check local mOTUs version) read profiling (Default: 0)
	-profileKraken [0|1]			1: perform kraken2 read proifiling (Default: 0)
	-estGenoSize [0|1]				1: estimate average size of genomes in data (Default: 0, currently not working)
	-krakenDB $PATH				$PATH to kraken2 database(s)

# IO for specific uses
	-newFileStructure [??]			just relink raw files for use in mocat
	-upload2EBI [??$PATH]				copy human read removed raw files to this dir, named after sample

# MODE: map2tar (map2DB / map2GC) mapping raw reads to reference databases (like genomes, functional DBs etc`)
# this mode is activated by calling ./MG-TK.pl map2tar -ref somthing.fa [..]
	-ref				reference database (.fa format)
	-mapUnmapped [0|1]			1: map unmapped reads (-saveReadsNotMap2Assembly) onto reference database
	-decoyMapping [0|1]				1: "Decoy mapping": map against reference genome AND against assembly of metagenome (drawing obvious better hits to metagenome, the "decoy") (Default: 1)
	-competitive2ndmap [-1|0|1|2]			1: Competitive, 2: combined but report separately per input genome, -1: combined and report all together. (Default: 1)
	-mapnms					name for files
	-redo2ndmap				
#D2s distance
	-calcInterMGdistance [0|1]		calculate nt distances between MGs; deprecated (no longer supported)
#Institute specific: EI
	-wcKeyJobs [#]				#: attach key to each job to help institute track cluster usage (defunct)
```

Comment: usually ‘0’ means switching a mode off, and ‘1’ means switching a mode on (unless specified).

</details>

### Genecat.pl flags

<details>
  <summary>Expand section</summary>

```{sh}
#Directories/files
	"GCd=s"  			Main save location for gene catalog and supporting files
	"tmp=s"				Tmp dir, global availalbe
	"glbTmp=s" 			Global tmp dir, same as -tmp usually
	"map=s" 			Mapping file, can be a combination of several .map files to combine different datasets (e.g. -map file1.map,f2.map)

#run modes
	"m|mode=s" 			possible modes: mergeCLs CANOPY specI kraken kaiju FMG_extr FOAM ABR FuncAssign protExtract ntMatchGC geneCat

#cluster options
	"clusterID=i" 			identity at which to cluster gene catalog, default: 0.95
	"minGeneL=i" 			minimal gene length for gene to be included in gene catalog, default: 100
	"extraGenesNT=s" 		add genes (nt) from external sources, e.g. from complete genomes
	"extraGenesAA=s" 		add genes (AA) from external sources, e.g. from complete genomes
	"mmseqC=i" 				1: use mmseqs2 instead of CD-HIT for gene clustering
	"decluterMatrix=i" 		1: declutering of gene matrix. Can give an edge to canopy based MGS, but also introduce unwanted biases. Default: 0

#flow control
	"1stepClust=i" 			Cluster incomplete genes separately? Default: 0
	"submitLocal=i"			Important run mode switch, to submit jobs while geneCat is runnning single core 
	"submSystem=s"			SGE, slurm submission systems
	"continue=i" 			Flow control, 1: continue with found files 0: delete existing (partial) gene cat, start again
	"cores=i" 			Number of cores being used
	"cores0=i" 			Specifcally cores only for the big main mmSeqs2 clustering job.. (takes a lot of mem and cores usually)
	"cores3=i" 			Num cores for small jobs that really don't require that much power.. 
	"mem=i" 				Max mem
	"mem3=i" => 			Max mem for smaller jobs
	"oldStyleFolders=i"		Deprecated. only used for results calculated with an older MG-TK version
	"sampleBatches=i"		How many batches to use for initial accumulation of genes? (200-500 samples per batch recommended). Default: Auto

#Binning/MGS related
	"binSpeciesMG=i" 		Use MAGs to create MGS? 1= metaBat2, 2=SemiBin, 3=metaDecoder
	"useCheckM2=i" 			1: use checkM2 completeness predictions, Default: 1
	"useCheckM1=i" 			1: use checkM completeness predictions, Default: 0
	"doStrains=i" 			1: calculate intraSpecific phylogenies on each MGS
	"doMags=i" 				1: start canopy clustering, metabat2 & subsequent merging into MGS
	"canopyAutoCorr=f" 		canopy clustering parameter to filter autocorrelated genes prior to canopy clustering

#Marker Genes/ taxonomy
	"MGset=s" 			Use either FMG or GTDB marker genes to compare and merge MAGs and calculate their abundance

#flags for specific modes
	"out=s" 			Output dir, only used in modes protExtract ntMatchGC 
	"functDB=s" 		for FuncAssign mode: functional DBs to annotate gene cat to 
	"refDB=s" 			For ntMatchGC mode: reference fasta DB 
	"fastaSplit=i"		For FuncAssign mode: split geneCat into chunks to parallelize jobs. Default: 500M 
```

</details>

### MGS.pl flags

<details>
  <summary>Expand section</summary>

```{sh}
	"GCd=s" 							#gene catalog dir
	"tmp=s" 							#temp dir
	"submit=i" 							#1:submit jobs, 0: dry run. Default: 1
	"canopies=s" 						#location of canopy clustering output file (clusters.txt)
	"smallCores=i" 						#cores used for normal jobs (not intensive)
	"bottleneckCores=i" 				#cores for compute intensive jobs
	"useRHClust=i" 						#1: do hierachical clustering of MGS genes. Default:0
	"redoRhcl=i"						#rewrite R hierachical clusterings
	"redoDeepCan=i" 					#rewrite deep corraltions to Rhcl clusters
	"redoTax=i" 						#rewrite tax annotations
	"MGset=s" 							#GTDB or FMG, which marker genes are used? Default: GTDB
	"mem=i" 							#memory used for intensive jobs
	"completeness=i" 					#what qual should final MGS have at least??
	"contamination=i"					#contamination threshold for accepting MGS
	"strains=i"							#1: calc instra species strain phylogenies. Default: 0
	"useCheckM2=i"						#CheckM2 default qual checking of MAGs/MGS
	"useCheckM1=i"						#CheckM default qual checking of MAGs/MGS
	"binSpeciesMG=i"					#0=no, 1=metaBat2, 2=SemiBin, 3: MetaDecoder
	"ignoreIncompleteMAGs=i" 			#1: assemblies without MAG calculations are ignored. Default: 1
	"legacy=i"							#1: use legacy code as pre Dec `22 (clustering is a bit more muddy, reported abundances slightly different, remember to use -MGset FMG). No longer supported. Default: 0
```
  
</details>

### buildTree5.pl flags

<details>
  <summary>Expand section</summary>

```{sh}
    #basic options
    -fna $PATH				path to .fna files
    -aa $PATH				path to .fa files
    -cats $PATH				path to category file, that sorts fna/aa sequences by a) functional category and b) species they originate from
    -map $PATH				map file
    -outgroup $PATH				path file with outgroup in sequence set
    -subsetSmpls [-1|#]				 default -1

     #options for phylogeny
    -runRAxML [0|1]				do RAXML (default 0)
    -runRaxMLng [0|1]				doRAXMLng (default 0)
    -runFastTree [0|1]				doFastTree (default 0)
    -AutoModel [0|1]				treeAutoModel, 1:iqtree: choose model automatically (a bit slower); (default 1)
    -SynTree [0|1]				default 0
    -NonSynTree [0|1]				default 0
    -bootstrap [0|1]				default 0
    -superTree [0|1]				doSuperTree (default 0)
    -superCheck [0|1]				doSuperCheck (default 0)

     #MSA related options
    -MSAprogram [0|1|2|3|4]				do MSA with 1:clustal or 0:msaprobs, 2:mafft (default), 3:guidance2, 4:MUSCLE5
    -minOverlapMSA [#]				min overlap in MSA columns, in order to retain column (default 0)
    -maxGapPerCol [#]				same as minOverlapMSA, but for MSAfix and %of gaps allowed in a column (default 1)
    -minPcId [0|1]				sequence is filtered from data, unless the average minPcId is >= minPcId; (default 0)

     #options for flow control
    -fixHeaders [0|1]				fix the fasta headers, if too long or containing not allowed symbols (nwk reserved) (default 0)
    -useEte [0|1]				using Ete (default 0)
    -calcDistMat [0|1]				distance matrix of either AA or NT (depending on MSA; default 0)
    -calcDistMatExt [0|1]				distmat of other AA or NT (depending on MSA), e.g. running two times an MSA; (default 0)
    -calcDiffDNA [0|1]				 default 0
    -postFilter $PATH				"," sep list of zorro,guidance2,macse

     #File operations
    -isAligned [0|1]				if input is already alligned (default 0)
    -rmMSA [0|1]				remove MSA, to save diskspace (default 0)
    -gzInput [0|1]				to save diskspace (default 0)

     #specific gene quality filters, useful for metagenomic data that is often incomplete
     -NTfilt [#]				fractions of nucleotides (NT) that need to be present in sequence to be included in final, combined MSA (default 0.8)
    -NTfiltPerGene [#]				if several genes represent a tree tip, the fraction of NTs that need to be present for the gene to be accepted in the final, combined MSA (default 0.1)
    -GenesPerSpecies [#]				min fraction of genes present after filtering, if below the species will be excluded from phylogeny (default 0.1)
    -fracMaxGenes90pct [#]				gene cats to keep, e.g. 25% of 90th percentile (default 0.25)
    -NTfiltCount [0|1]					total NT count (default 0)
    -runLengthCheck [0|1]				check that sequence length can be divided by 3

     #popgen related options
    -runClonalFrameML [0|1]				doCFML
    -runGubbins [0|1]				doGubbins (default 0)
    -runDNDS [0|1]				run dNdS analysis
    -runTheta [0|1]				doTheta
    -genesForDNDS $PATH				file with list s{,} with selected genes just for dnds
    -DNDSonSubset [0|1]				run dnds just on subset (given by genesForDNDS) of genes
    -codemlRepeats [#]				repeatCounts, set how often each model should be repeated to check for convergence (default 2)
    -outDCodeml $PATH				codemlOutD,
    -genesToPhylip [0|1]				doGenesToPh (default 0)
    -runFastgear [0|1]				doFastGear (default 0)
    -runFastGearPostProcessing [0|1]				doFastGearSummary (default 0)
    -clustername $PATH				clusterName
```
  
Comment: usually ‘0’ means switching a mode off, and ‘1’ means switching a mode on (unless specified).
  
</details>











## Outputs

### General output structure

In MG-TK, output can be either in 

a) the output dir defined in the .map file. These are usually 1) assemblies, and MAGs, gene predictions, abundances of these 2) SNP callings (==corrected assemblies for practical purposes) 3) assembly independent methods such as riboFinder, metaPhlan, mOTUs, mappings of reads against assemblies/reference databases, 4) outputs from map2tar mode.

b) gene catalogue related outputs are stored exclusively in the output folder defined in your geneCatalog script

c) MGS (metagenomic species) related outputs, representing a merge of gene catalogue and assemblies/MAGs, are by default stored in the  `gene_catalog_ouput_folder/Bin_SB/`. However, output can be changed with the `MGS.pl -outD X` flag.

### Abundance Matrices

Several feature abundance tables are produced for both functional and taxonomic abundances. Tables are all in tab-separated format, oriented with features on rows and samples on columns. For features which which have a hierarchical structure, matrices summing the features at each level of the hierarchy are provided and usually denoted in the file names as {..}L0.txt, {..}L1.txt through the different hierachical levels.

#### Taxonomic

Feature names are the GTDB lineage up to a given rank, separated by semicolons. When taxonomy could not be identified at a certain rank a ? is used. For MGS level taxonomy, the identifier of the MGS is appended to the full lineage. A feature `-1` or `?` is included often as the first row a matrix which represents genes which could be assigned to any taxonomy. `-1` usually represents reads that could not be matched at all, while `?` represents hits where the taxonomy/functional assignment is unclear for various reasons.
 - `Bin_SB/Annotation/Abundance`: Abundance based on clustering of the gene catalog into MGS. Abundance is estimated by the coverage of *core genes* included in the MGS. Matrices provided at ranks from domain (`MGS.matL0.txt`) to species (`MGS.matL6.txt`), as well as individual MGS (`MGS.matL7.txt`). At MGS level, some features end with ? rather than and MGS identifier, these are taxa for which no MGS was present; they are identified solely through the LCA assignments of marker genes.
 - `Anno/Tax/GTDBmg_MGS`: **Deprecated, do not use any longer!!** Abundance based on both clustering of the gene catalog into species provided by GTDB (or specI, from proGenomes, if the `FMG` option is used). Abundance is estimated by the coverage of *all genes* includes in the MGS. Matrices are provided at ranks from superkingdom (`specI.superkingdom`) to species (`specI.species`), as well as at individual MGS (`specI.mat`). In the MGS matrix, the entries which represent MGS do not have their full lineage output; the lineage can be found in the file `MGS2speci.txt` or in `specI.tax`, the latter of which includes lineages for non-MGS entries in the `specI.mat` matrix.
 - `Anno/Tax/GTDBmg`: **Deprecated, do not use any longer!!** As above, but based only on specI and not on MGS clustering. In this case, the `specI.mat` and `specI.species` matrices are identical.

 #### Functional
 Abundance of functional annotations using several databases are available in `Anno/Func`. Abundance is based on coverage of genes with a given functional annotation regardless of taxonomy. For those with hierarchical structures, data is provided either summed to each level of this hierarchy, or for the lowest level of hierarchy and the structure provided in a separate file, as noted below.

**KEGG, SEED, BSB modules**
 The modules are in `Anno/Func/modules` (KEGG modules in the directory `/modules`). In each of these directories are:
 - `KEGG.mat`: Abundance matrix, giving abundance of each module in each sample.
 - `KEGG.descr`: Table describing the hierarchy. Module identifier in the first column, then name at each level in following columns.
 - `KEGG.KOused`: The individual functions identifiers (i.e. KEGG orthologs) which were members of this module in each sample.
 - `KEGG.MODscore`: Completeness of the module for each sample (proportion of members which are observed at least once).
 
 The abundance of KEGG Orthologs is provided in `Anno/Func/KGML0.txt`
 
 **CAZy, eggNOG, TCDB**
 For these databases, abundance matrices are provided summed to each level in `Anno/Func`, with a prefix indicating the database, and `L{i}.txt` the level, with 0 being the most specific, e.g. `CZyL0.txt` giving the most specific level of CAZy annotation abundance.
 
 **gene assignments**
It can be useful to know which functional annotations a given gene has, or vice versa. This can be found in `Anno/Func/DIAass_{db}.srt.gzgeneAss.gz`, where `{db}` is the databse of interest (KEGG, CAZy etc.) This file gives the gene id in the gene catalogue as the first column, followed by it's annotation at the lowest level of the database. A gene can have more than one annotation in this file.

### MAG/MGS Gene Content
The gene content of a MAG can be found in `Bin_SB/LOGandSUB/MAGvsGC.txt.gz`. The columns *MAG* and *MGS* give the MAG id and the MGS which it is part of. The column *Representative4MGS* contains a star when the MAG is the representative for this MGS. Each marker gene has a column giving the genes from the MAG which matched against that marker. 

The final column, *other_genes*, gives a comma separated list of all the other genes which are part of this MAG. These are in the order they appear on the contig, with a double comma indicating the start of a new conting. Using this list and the gene function table explained above lets you link MAGs/MGS and function.

