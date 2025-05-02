
## Running MG-TK

MG-TK is programmed for HPC environments (Linux) and was conceptualized to process 1000's of metagenomes. It relies therefore on job schedulers (slurm, SGE and LSF are supported) and multiple safeguards to resume failed jobs. Please see examples below for specific runs.

### Temporary and output files

The output path for storing non-temporary files (like assemblies, binnings, gene predictions etc), is defined in each mapping file separately, composed of the arguments "#OutPath" and "#RunID". The final output will be stored in the dir "#OutPath/#RunID/", here each sample will have its own folder, and within this folder assemblies, gene predictions (assembly dir), mapping reads to the assemblies (mapping dir) and a detailed log of the steps run (LOGandSUB dir), will be stored.

Since the pipeline is expected to run on a compute cluster, temporary directories are of enormous importance for a) performance and b) file exchange between compute nodes that are usually physically separated clusters.
The pipeline expects a path to a storage that is globally available on all nodes and a tmp dir that is locally available on each node (given by arguments "globalTmpDir" and "nodeTmpDir" in the config file). 

### Mapping file for MG-TK


Most importantly you need a mapping file to your files. See 'examples' dir for some map examples (also explaining how to do compound assemblies, compound mapping). These column names (headers) are reserved key words in the mapping file (other columns can be eg. metadata per sample etc):
- **#SmplID** [STRING] MG-TK maps always need to have the first column names *#SmplID*. The string in this column will used in all subsequent analysis, intermediary files, sequence heads etc to uniquely identify samples, therefore choose with extreme care! Good practice would be to include some basic information about the sample in the SMPLID, but should be as short and descriptive as possible. *DO NOT USE SPECIAL CHARACTERS IN THE SMPLID, keep it basic*!  
- **Path** [STRING] - is the relative path to fastq[.gz] files for each sample (see #DirPath, this needs to be set to the absolute path). All files ending with .fq or .fastq (can have .gz after) in the dir will be used for that specific samples. 1. or 2. indicates first or second read. E.g. al0-0_12s005629-2-1_lane3.2.fq.gz is the second read, here the pipeline expects to have al0-0_12s005629-2-1_lane3.1.fq.gz in the same dir.  
Further, you can add the following specifics for each single sample:   
- **AssmblGrps** [STRING] - set this to a number or string. all samples with the same tag will be assembled together (e.g. samples from the same patient at different time points).  
- **MapGrps** [STRING] - set a tag here as in AssmblGrps. All reads from these samples will be thrown together, when mapping against target sequences (only works with option "map2tar" and "map2DB").
- **SupportReads** ['PB', 'mate'] - in case you have additional reads, that are not normal illumina hiSeq, e.g. miSeq or hiSeq in mate pair sequence mode ('mate') or PacBio reads ('PB').
- **SeqTech** ['ill', 'ONT', 'PB', 'SLR'] - Sequencing technology used in sample: illumina short reads ('ill') , Oxford Nanopore ('ONT'), PacBio ('PB') or sythetic long reads ('SLR').
- **ReadLength** - Expected read length in sample. Is usually automatically determined, use with caution!
- **EstCoverage** [0/1] - (Deprecated!!) Used to indicate if the avg coverage of genomes should be estimated in sample.
- **SupportReads** [tag:path] - Additional reads created with a different seq technology. E.g. miSeq ('miSeq:/path/to/file'), mate-pair ('mate:/path/to/file') or PacBio ('PB:/path/to/bam').
- **ExcludeAssembly** [0/1] - Exclude sample from assemblies?
- **cut5PR1** [INT] - remove the first nts (from 5') on read 1
- **cut5PR2** [INT] - remove the first nts (from 5') on read 2
- **firstXreadsRd** [INT] - only read the first X reads (paired reads count as 2) for that sample.
- **firstXreadsWr** [INT] - only write the first X reads (paired reads count as 2) for that sample.
 
 The following tags can be added to a new line (ie row) in the map. Tag is followed by tab delimiter and specific input.

#### Required map tags
- **#OutPath**	[Path] Where to write the output (can be massive, make sure you have enough space)
- **#RunID**	[string] The directory below OutPath, where results are stored. Also serves as global identifier for this run
- **#DirPath**	[Path] Base directory where subdir with the fastqs can be found. You can insert this on several lines, if the base path changes for all samples afterwards.

#### Optional map tags

- **#NodeTmpDir**	[Path] temporary dir only accessible within each compute cluster node, overrides **nodeTmpDir** definition in config file
- **#GlobalTmpDir**	[Path] temporary dir (scratch) accessible from all compute nodes, overrides **globalTmpDir** definition in config file
- **#mocatFiltPath**	If for some reason you are forced to use mocat filtered fastqs and not the original, unfiltered files (strongly recommended), than you can indicate in which subdir these mocat files can be found
- **#RelaxSMPLID**	[TRUE/FALSE] 	Use FALSE to deactivate basic checks if the #SmplID adheres to MG-TK formats. Caution: use on your own risk!
- **#WARNING**	[OFF/ON]	If **OFF** MF won't stop when an error is encountered in the map. Caution: use on your own risk!

After this follow the sample IDs and the relative path, where to find the input fastqs.  
See _examples/example_map_assemblies.map_ for a very complicated mapping file with several source dirs.

#### Example mapping file

```{sh}
#SmplID	Path	SmplPrefix	AssmblGrps
#OutPath	/hpc-home/path/to/your/results/folder
#RunID	NAMEofresultsFOLDER (#MG-TK will make this folder with this name by itself)
#DESCRIPTION (#not important, but you can mark what is the run is about)
#DirPath	/path/to/folder/with/raw/reads
Mouse11t0		PID_C11T0_	M11
Mouse11t1		PID_C11T1_	M11
Mouse12t0		PID_C12T0_	M12
Mouse12t1		PID_C12T1_	M12
Mouse14t0		PID_C14T0_	M14
Mouse14t1		PID_C14T1_	M14
Mouse15t0		PID_PD11T0_	M15
#DirPath	/path/to/another/folder/with/more/raw/reads
Mouse15t1	SubDir1		M15
Mouse16t0	SubDir2		M16
Mouse16t1	SubDir3		M16
```		

#### Tips and recommendataions for creating mapping files

- It is recommended to create the mapping file in **Excel** and copy-paste it in a **.map** text file afterwards (will be tab-delimited by default, the expected MG-TK format). You can use functions like "=VLOOKUP()" to match sample IDs across different tables. 

- The **#SmplID** column determines the name of a sample all the way throught the pipeline! Be very careful what ID you choose, as this will impact the sample names you'll have to deal with later, choose something a) short and b) descriptive. Avoid c) special characters (_|$%~\`\*& etc) in the SmplID!

- **AssmblGrps**: Assembly groups are useful for assembling samples from e.g. a time series together, giving a better assembly usually. Choose the name of an assembly group a) unique b) short and descriptive and c) avoid special chars (\[\]{}_|$%~\`\*& etc)!

- If using **assembly groups**, try to keep samples from the same assembly group as a block. MG-TK can also deal with these assembly groups distributed across the map, but in terms of job submission strategy it's best to have these samples next to each other in the map (and also for you organizing your experiment).

- <ins>**Loading and saving a mapping file into R will likely lead to problems!**</ins> This is because the #DirPath tag sets the path for all samples underneath. Loading this into R will often skip the #DirPath line or reorder the samples, so saving this again will lead to wrong paths being set!







## Additional usage scenarios


MG-TK can be used for a bulk of tasks not directly related to initial assembly, profiling or MAGs, but often related to postprocessing these. Two usage scenarios (map2tar and building phylogenies) are listed below.


### map2tar mode

- This mode maps the reads to `reference` genomes, e.g. from a mock community. This allows a fast profiling for specific purposes. The mode switches off assembly-based processes.

1. create mapping file `path/to/mapping_file.map`

```{sh}
#SmplID	SmplPrefix	AssmblGrps
#OutPath	/path/to/results/profilingMF
#RunID	name_of_the_run_dir
#DESCRIPTION
#DirPath	/path/to/raw/reads
BERG100	BERG100	BERG100
BERG100w	BERG100w	BERG100w
BERG10	BERG10	BERG10
BERG10w	BERG10w	BERG10w
BERGmock	BERGmock	BERGmock
```

2. make script with first command `run_profiling.sh`

```{sh}
MAP="/path/to/mapping_file.map"
perl $MF3DIR/MG-TK.pl map2tar \
	-map $MAP  -ref '/path/to/reference/mock_community/*.fasta' -filterHumanRds 0 -mappingCores 12 -mapperFilterIll '0.02 0.75 00'  -redo2ndmap 0 -mappingMem 15 -submit 1 -competitive2ndmap -1 -decoyMapping 0
```

- Explanation: ref are the .fasta formated reference genomes you want to map your metagenomic reads to, metagenomic reads are defined in the map, as in other runs. -mapperFilterIll defines how the mapped reads will be quality filtered. -competitive2ndmap defines if reads will be mapped against all references at once (competetitive) or separately against each single reference. -decoyMapping determines if an already created read assembly will be used to "decay" map reads against (useful if you suspect that most reads aligning to your reference would be false positive hits).

- run with `bash run_profiling.sh`

- output files will contain coverage per window, contig, etc. which can be used for plotting.

### Building phylogenetic trees with MG-TK

1. create a file `phyloScript.pl`

```{sh}
#!/usr/bin/perl
use strict; use warnings;

my $bts = "/path/to/MG-TK/secScripts/phylo/buildTree5.pl";
my $inD = "/path/to/input/dir/phylo/";
my $outD = $inD."/bts/";
my $tempD = "/path/to/scratch/dir/treetest/";
my $numCore = 8;

my $cmd = "$bts -genoInD $inD -outD  $outD -tmpD $tempD -runIQtree 1 -iqFast 1 -AAtree 1 -cores $numCore -wildcardflag '/*.f*a' -continue 1 \n";

print $cmd;
system $cmd;
```

Explanation: $inD is an input dir with complete genomes, the script will extract FGMs and build tree between genomes. `-AAtree 1` tells the script to use AA MSAs to build the phylogeny via iqTree. `-wildcardflag '/*.f*a'` tells the script how to look for reference genomes in $inD. 

2. run the script `phyloScript.pl`

- Run with `perl /path/to/phyloScript.pl` together with submission script on cluster.

3. you can do many additional phylogeny / pogen related analysis with the `buildTree5.pl` script ([see flags](#buildtree5.pl-flags)) 

