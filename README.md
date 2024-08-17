# MG-TK

## Table of Contents

<details>
  <summary>Expand Table of Contents</summary>
  
- [Introduction](#introduction)
- [Requirements](#requirements)
- [Installing MG-TK](#installing-MG-TK)
- [Running MG-TK](#running-MG-TK)
	- [Temporary and output files](#temporary-and-output-files)
	- [Mapping file](#mapping-file)
- [Examples](#examples)
	- [MG-TK metagenomic assembly and gene catalog](#MG-TK-metagenomic-assembly-and-gene-catalog)
	- [Assembly-independent MG-TK mode](#assembly-independent-MG-TK-mode)
	- [Hybrid assemblies](#hybrid-assemblies)
- [Outputs](#outputs)
	- [Abundance matrices](#abundance-matrices)
	- [Gene function & MAG/MGS gene content](#gene-function--magmgs-gene-content)
- [Additional information](#additional-information)
	- [MG-TK.pl flags](#MG-TKpl-flags)
	- [Genecat.pl flags](#genecatpl-flags)
	- [MGS.pl flags](#mgspl-flags)
	- [buildTree5.pl flags](#buildtree5pl-flags)
- [Additional usage scenarios](#additional-usage-scenarios)
	- [map2tar mode](#map2tar-mode)
	- [Building phylogenetic trees with MG-TK](#building-phylogenetic-trees-with-MG-TK)
- [Trouble shooting](#trouble-shooting)
- [License, citations etc](#license,-citations-etc)

</details>

## Introduction 

MG-TK is a pipeline developed to 
- Assemble metagenomes, profile miTags, profile functions, profile taxonomy, build MAGs (metagenomic asssembled genomes) using a variety of approaches (MG-TK.pl)
- Build a gene catalog based on these assemblies and predicted genes, build abundance matrices from these and annotate the genes functionally (secScripts/geneCat.pl)
- Merge MAGs into MGS (metagenomic species), build inter- and intra-speciies phylogenies (secScripts/MGS.pl)
- Additional functionalities are available, such as building phylogenetic tree automatically for many genomes and calcualting population genetic statistics on these (secScripts/phylo/buildTree5.pl)

MG-TK is implemented in Perl, C++ and uses R and python scripts. 

Author: Falk Hildebrand <Falk.Hildebrand@gmail.com>
## Installing MG-TK

<details>
  <summary> Expand Installation section </summary>
  

### Requirements

MG-TK requires a perl installation and sdm requires a fairly recent C++ compiler (like gcc or clang) that supports C++11.
MG-TK currently only works under linux, and is expected to run on a computer cluster. Since the pipeline includes a lot of external sofware, you will need fully installed Micromamba ([https://mamba.readthedocs.io/en/latest/installation.html](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)).


### Installation


MG-TK can be downloaded directly from Github, using:
```
git clone https://github.com/hildebra/mg-tk.git
```
MG-TK comes with an installation script, that uses micromamba. Ensure you have micromamba installed for your account on a linux HPC. Then run:
```
bash helpers/install/installer.sh
```

This will guide you through the installation (should run completely automatic) and requires internet access. Since a lot of packages will be installed, this can take an hour or longer. All required software will be downloaded and installed in the Conda/Mamba directories.

If you are having issues with package conflicts when `installer.sh` is creating environments, trying setting your channel priority to flexible: `micromamba config set channel_priority flexible`

Last, you can run 
```
./MG-TK.pl -checkInstall
```
to check that some essential programs have been correctly installed and are available in the exptected environments. Note that this is only a subset of programs, but should cover most use cases of MG-TK. (This will also automatically run after each installation of MG-TK)


### Updating MG-TK

MG-TK will be frequently updated. To get the latest version, go to your MG-TK directory and run
```
git pull
```
Sometimes new packages will be included or program versions modified. To obtain these changes, run the install script again (this will update existing environments - no worries, this is not a complete reinstall):
```
bash helpers/install/installer.sh
```

### Preparing MG-TK

- follow installation process (essentially `git clone https://github.com/hildebra/mg-tk.git` & run `bash helpers/install/installer.sh` )

- After the instalation is complete, you will find the file named: "config.txt" inside of the MATAF3 directory. This is the main file where you have defined all the paths for directories and slurm configuation. Always check in order to ensure that all directories are correct: 

    - MFLRDir	`/path/to/your/mg-tk/installation/`
    - DBDir	`/path/to/your/database_dir/`

- change tmp dir (scratch space) with project scratch folder:

    - globalTmpDir	`/path/to/your/scratch/`
	- nodeTmpDir	`/path/on/node/to/tmp` -> on slurm systems this could be a variable, e.g. `$SLURM_LOCAL_SCRATCH/MG-TK/`

- follow either assembly-dependent or assembly-independent tutorial

</details>


## Running MG-TK

MG-TK is programmed for HPC environments (Linux) and was conceptualized to process 1000's of metagenomes. It relies therefore on job schedulers (slurm, SGE and LSF are supported) and multiple safeguards to resume failed jobs. Please see examples below for specific runs.

### Temporary and output files

The output path for storing non-temporary files (like assemblies, binnings, gene predictions etc), is defined in each mapping file separately, composed of the arguments "#OutPath" and "#RunID". The final output will be stored in the dir "#OutPath/#RunID/", here each sample will have its own folder, and within this folder assemblies, gene predictions (assembly dir), mapping reads to the assemblies (mapping dir) and a detailed log of the steps run (LOGandSUB dir), will be stored.

Since the pipeline is expected to run on a compute cluster, temporary directories are of enormous importance for a) performance and b) file exchange between compute nodes that are usually physically separated clusters.
The pipeline expects a path to a storage that is globally available on all nodes and a tmp dir that is locally available on each node (given by arguments "globalTmpDir" and "nodeTmpDir" in the config file). 

### Mapping file

<details>
  <summary>Expand section</summary>

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


- Naming: Use simple names as SmplID that are self descriptive. Avoid special characters, that is try to stay with letters in \[a-zA-Z0-9.\] character space for the smpl ids

</details>
 

## Examples

<details>
  <summary> Expand example assembly with MG-TK </summary>
  
### MG-TK metagenomic assembly and gene catalog

The figure below shows example of steps involved in the assembly-dependent mode. White rectangles indicate inputs and outputs, grey boxes name each of the steps, and yellow boxes show names of the scripts that are generated and submitted in each step. Blue boxes indicate additional steps that are required for subsequent MGS analysis.

<img src="./helpers/documentation/assembly-dependent.svg" style="width: 800px;"/>


#### 1. create mapping file

Typically you would use Excel to create the mapping file and copy-paste it later into a text file (will be by default tab-delimited). This text file, typically with the file ending **.map** can then be saved to the HPC.

#### 2. make script with first command `RUN.mfc`

Insert your MG-TK command, a bash and slurm header in RUN.mfc. 

Example:

```{sh}
#!/bin/bash
#SBATCH -J SUB_MF
#SBATCH -N 1 --cpus-per-task=1 --mem=10024 --export=ALL
#SBATCH -o [currentDir]/run_mataf3_mhit.mfc.otxt
#SBATCH -e [currentDir]/run_mataf3_mhit.mfc.etxt
#SBATCH -p "ei-long,qib-long"

set -e
ulimit -c 0;
MAP=/path/to/your/mapping/file/FILE.map
perl $MF3DIR/MG-TK.pl -map $MAP-assembleMG 2 -spadesCores 12 -spadesKmers "25,43,67,87,111,131" -spadesMemory 100 -mapReadsOntoAssembly 1 -kmerPerGene 0 -filterHostRds 1  -filterHostKrak2DB /hpc-home//data/DB/kraken2/hsap/ -mappingMem 5 -profileMOTU2 0  -profileMetaphlan3 1 -Binner 2  -maxConcurrentJobs 600 \
-from 0 -to 1 -submit 1 -getAssemblConsSNP 0
```

- MG-TK now does

    - read filtering and cleaning
    - read profiling (if set like in the above command)
    - assembly per assembly group
    - mapping to the assemblies


#### 3. Test the run and the map:

- always run with submit 0 before you finally decide on parameters and final run

- -from X -to Y controls that only samples X to Y will be processed. Good for testing e.g. only first sample in map (-from 0 -to 1 ). Setting -to to very high number will just run to the end of the map (e.g. -to 99999 to finish map).

- if you want to run strain analyses later, set `-getAssemblConsSNP 1` to calculate consensus SNPs needed for strain analysis

#### 4. Running MG-TK

- Run with `bash run_mataf3_mhit.sh` --> this will submit a lot of different jobs to the HPC queue

 console output:
```
        This is MG-TK 0.33
        Using qsubsystem: slurm
        Using qsubsystem: slurm
        /projects/data/results/mataf3_test1/LOGandSUB/qsub.log
        Reset range of samples to 40

        ======= Mouse11T0 - 0 - M11T0 =======
        1:2  1:0
        SUB:_UZ0        SUB:_SDM0       SUB:_cln0       Running Contig Stats on assembly

        ======= Mouse11T1 - 1 - M11T1 =======
        2:2  1:0
        SUB:_UZ1        SUB:_SDM1       Assembly stepSUB:_A1    SUB:_GP1        SUB:_cln1       Running Contig Stats on assembly
        SUB:_CS1
```

- more advanced usage: run `sbatch run_mataf3_mhit.sh`. This will submit the job to the cluster queue, and from there the .mfc job will submit more jobs. The output from MG-TK will be stored in `#SBATCH -o [currentDir]/run_mataf3_mhit.mfc.otxt` and `#SBATCH -e [currentDir]/run_mataf3_mhit.mfc.etxt` defined above.


#### 5. rerun MG-TK

MG-TK is conceptualized to detect automatically if certain steps need to be run again (e.g. because the job crashed or some files from other subjobs were not yet available). Therefore you will usually need to rerun the same MG-TK command several times. **However, before restarting MG-TK make sure that all job submissions from your previous run have completed (or don't start due to job dependecies)!**

- Once MG-TK detects no further jobs to be submitted, it will let you know (check the output of the RUN.mfc command). At this point you can advance to creating a gene catalog. MG-TK will create a `GeneCat.sh` script (detailed in RUN.mfc output). 

#### 6. building a gene catalog by running the `GeneCat.sh` script

- After every sample has successfully passed through the pipeline, MG-TK produces GeneCat_pre.sh script that needs to be adapted:

- In the `GeneCat.sh` script you need to specify an output directory, max memory usage and the cores you want to use. After that the script can look like this:

```{sh}
#!/bin/bash
#SBATCH -N 1 --cpus-per-task=1
#SBATCH -o /ei/projects/data/results/GeneCat_pre.sh.otxt -e /ei/projects/data/results/GeneCat_pre.sh.etxt
#SBATCH --export=ALL --mem=81920 -J myFirstGeCat
#SBATCH -p "ei-medium,qib-medium,ei-long,qib-long"
set -e
ulimit -c 0;

#creates gene catalog in the specified outdir with specified cores, attempting to reuse existing dirs (in case catalog creation failed):
perl /hpc-home/project/MATAF3/secScripts/geneCat.pl \
		-map /ei/projects/data/results/mapping_file.map \
		-GCd /ei/projects/data/results/genecat \
		-mem 200 -cores 24 -clusterID 95 -doStrains 0 -continue 1 \
		-Binner 2 -useCheckM1 0 -useCheckM2 1 -MGset GTDB 
```

- now you have to run the GeneCat.sh script with `sbatch GeneCat.sh`, which will submit a lot of different jobs to the cluster.

- This step of the pipeline calls `GeneCat.pl` and does:

    - creates a gene catalog (at 95% id) from predicted genes using mmSeqs2
    - extracts proteins corresponding to genes in gen catalog
	- identifies genes that represent marker genes (GTDB)
	- creates a gene abundance matrix (literally gene catalog millions of genes and their abundance in all metagenomic samples.. very big files!)
    - assigns basic functions to genes in gene catalog (KEGG, eggNOG, CAZy, ARG, ..)
    - accumulates MAGs and dereplicates these into MGS (metagenomic species)
	- calculates abundances of MGS in samples and an MGS taxonomy
	- calculates intra-specific phylogenies for each MGS (basically strain tracking across metagenomes)

</details>


### Assembly-independent MG-TK mode

<details>
  <summary>Expand example assembly-independent MG-TK mode </summary>


This mode is especially helpful if your metagenome is from a highly complex (e.g. Soil) microbial community. In such cases often assembly is not possible, instead metagenomic reads can still be mapped to a reference.

#### 1. create mapping file:

	the same as for assembly mode `mapping_file.map`

#### 2. make script with first command `run_independent.mfc`

```{sh}
#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -o /hpc-home/project/run_independent.mfc.otxt
#SBATCH -e /hpc-home/project/run_independent.mfc.etxt
#SBATCH --mem=102400
#SBATCH --export=ALL
#SBATCH -p "ei-medium,ei-long,qib-medium"
#SBATCH --time=24:00:00
# #SBATCH --localscratch=ssd:10
#SBATCH -J SUB_MEST
set -e
ulimit -c 0;

#ABRc,CZy,PAB,NOG

MAP= /path/to/dir/mapping_file.map

perl $MF3DIR/MG-TK.pl -map $MAP -inputFQregex1 '.*_1\.f[^\.]*q\.gz$' -inputFQregex2 '.*_2\.f[^\.]*q\.gz$' \
	-mergeReads 0 \
	-profileFunct 1 -reParseFunct 0 -reProfileFunct 0 \
	-diamondDBs KGM,NOG,CZy \
	-diamondCores 8 \
	-maxConcurrentJobs 300 \
	-profileRibosome 1 -reProfileRibosome 0 \
	-profileMOTU2 0 -profileMetaphlan3 1 \
	-filterHostRds 0 -inputReadLength 150 -assembleMG 0 \
	-submit 1 \
	-from 0 -to 40
```

#### 3. run the script `run_independent.mfc`

- Run with `sbatch river_independent.mfc` --> this will submit a lot of different jobs to the HPC queue.

#### 4. rerun MG-TK  

- wait till all current jobs are finished, then rerun `sbatch run_independent.mfc`. This will check all jobs completed, and if so, create the feature abundance tables requested. In the above case these are KEGG, eggNOG and CAZy functional tables, as well as metaphlan3 and miTAG taxonomic tables.

</details>

### Hybrid assemblies
<details>
  <summary>Expand example Hybrid assemblies MG-TK mode </summary>

With hybrid assemblies we mean metagenomic samples that were sequenced with multiple sequencing technologies, such as ill+ONT, ill+mate pair ill or ill+PB. Currently MF support ill+ONT in a basic form (ONT reads only used to improve assemblies) and a more advanced form for ill+PB. This section is describing the ill+PB mode, as it is more powerful, but also more complicated to use.

#### 1. setup ill+PB hybrid mode

Create a .map as usual for the illumina reads. Afterwards, add the following column to the .map: **SupportReads**
For each samples, where support reads (i.e. the PacBio reads) exist, add them to the **SupportReads** column, by indicating that they are PacBio reads through the prefix **PB:** followed by the path to the **.bam** (PB reads are usually saved in bam format). If you have several .bams for the same sample, these could be comma separated and will be automatically used by MF.

E.g. the map could now look like:

```
#SmplID	SmplPrefix	SupportReads	INFO	SeqTech	Path	AssmblGrps	sampleID	Individual	TimePoint
#RunID	PB.PAGE2								
#WARNING	OFF
#OutPath	/hpc-home/hildebra/grp/data/projects/								
#DirPath	/ei/projects/8/88e80936-2a5d-4f4a-afab-6f74b374c765/data/cloudpool/data/raw/Public/PRJNA529586/
S4zm	SRR8797713	PB:/path/tp/PB//test.hifi_reads.bc1011_tmp.bam		hiSeq				S4zm_R1177-S0001	0
S4qia	SRR8797712	PB:/path/tp/PB//test.hifi_reads.bc1012_tmp.bam,/path/tp/PB//test.hifi_reads.bc1011_tmp.bam		hiSeq				S4zm_R1177-S0001	0
S4zmGG1	SRR8797713	PB:/path/tp/PB//test.hifi_reads.bc1011_tmp.bam		hiSeq		I1		S4zm_R1177-S0001	0
S4qiaGG2	SRR8797712	PB:/path/tp/PB//test.hifi_reads.bc1012_tmp.bam		hiSeq		I1		S4zm_R1177-S0001	0
S4qiaGG3	SRR8797712	PB:/path/tp/PB//test.hifi_reads.bc1012_tmp.bam		hiSeq		I1		S4zm_R1177-S0001	0
S4zmS1	SRR8797713			hiSeq		I2		S4zm_R1177-S0001	0
S4qiaS2	SRR8797712			hiSeq		I2		S4zm_R1177-S0001	0
S4qiaS3	SRR8797712			hiSeq		I3		S4zm_R1177-S0001	0

```
(note that this map shows a mix of assembl grps, of samples with and without support reads. Further note that the tag **#WARNING	OFF** has to be used, since in this test case samples are being reused - something that would normally trigger MF to stop the run.)

#### 2. setup MF run

Setup your MF run like you would normally setup an assembly dependent metagenomic analysis (see examples above). However, these flags should be defined:
 - `-mapSupportReadsOntoAssembly 1`: setting this to `1` will lead to support reads (PB reads in this case) being mapped onto the assembly. A coverage profile is created that is then **separately** used from the illumina coverage profile in the binning step (which in our experience can significantly boost the recovery of MAGs).
 - `-mapper -1`: set the choice of mapper to the default (`-1`) which means that a mapper will be automatically selected. You could set this to e.g. **1** to do all mapping with bowtie2, or **3** to use minimap2 everywhere, but **-1** is the recommended choice as this will use bowtie2 for ill reads and minimap2 for PB reads.
 - `-assembleMG 5`: this flag is crucial and tells MF to conduct a hybrid assembly, using **megahit for illumina** and **metaMDBG for PB** reads.
 - `-inputReadLengthSuppl 8000`: not crucial, but good to have. Here we estimated that our support PacBio reads are on average 8000 bp long.
 
Your example MF call could now look like:
 
```{sh}
MAP=/path/to/map/PB_hybrid.map

perl $MF3DIR/MG-TK.pl -map $MAP -inputFQregexSingle '.*\.fastq\.gz' -inputFQregex1 '(.*_R1_001\.fastq\.gz)|(.*[_\.]1\.f[^\.]*q\.gz)$' -inputFQregex2 '(.*_R2_001\.fastq\.gz)|(.*[_\.]2\.f[^\.]*q\.gz)$' -inputBAMregex '.*\.bam$' \
-assembleMG 5 -spadesCores 12 -spadesKmers "25,43,67,87,111,127" -spadesMemory 200 -MetaBat2 2 \
-mapper -1 -mapSupportReadsOntoAssembly 1 \
-filterHostKrak2DB /path/to/kraken2/hsap/ -filterHostRds 1 \
-getAssemblConsSNP 1 -rmSmplLocks 0 \
-submit 1 -inputReadLength 150 -inputReadLengthSuppl 8000 \
-from 0 -to 1
```

#### 3. running MF in hybrid mode

We opted for a bit of a complicated processing of hybrid assemblies, that in the end allows for both usage of existing paths in MF as well as supporting more complicated (assemblyGrps) sample setup.
What this means for you as the user is simply that you need to run the above MF command several times. In the first iteration this will trigger the megahit illumina-assembly, in the second iteration remaining samples in assemblyGrp are mapped onto this illumina-assembly, in the third iteration the hybrid assembly with metaMDBG is started, fourth iteration will map remaining samples in assemblyGrp onto hybrid-assembly, fifth iteration will then finally starting consensus SNP calling and binning. Remember to wait between iterations until all current jobs have finished (though sample locks should normally prevent double submissions).

_tldr; this mode requires several iterations to complete_ 


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



## Additional information

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
	-maxConcurrentJobs			max jobs in queue, useful for large samples sets, currently only works on slurm 
	-excludeNodes [string]					exclude certain HPC nodes, comma separated list e.g. node1,node2,..
	-submSystem [qsub,SGE,bsub,LSF]	set submission system (default: autodetect)
	-redoContigStats [0|1]				if any step of requested analysis failed, contigStats (coverage per gene, kmers, GC content) will be deleted & started again
	-loopTillComplete [X:Y]			script will loop over the assigned samples until all jobs are finished #use synatx "X:Y" where X is num loops, Y is the window size, eg "6:250" would run 6 loops of max 250 samples, then move on to next 250 samples (#dangerous flag)
	-requireInput [0/1]		in case input reads are not present (e.g. something wrong in map), 0 will continue pipeline, 1 will abort
	-silent [0/1] 			Controls how much information is printed on console
	-OKtoRWassGrps [0|1]			1: can delete assemblies, if suspects error in them, powerful, but careful! (Default: 0)

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
	-filterHumanRds [0|1]			0: do not filter host reads, 1: filter host reads (same as -filterHostRds)
	-filterHostRds				0: do not filter host reads, 1: filter host reads (same as -filterHumanRds)
	-filterHostKrak2DB $PATH		path to host kraken database
	-onlyFilterZip [0|1]			??
	-mocatFiltered [0|1]			??
	-filterHostKr2Conf [#]                      set host kraken2 confidence parameter (e.g. 0.05) 
	-filterHostKr2Quick ["--quick"]             set quick option for kraken2 (should be "--quick")

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

# binning
	-Binner [1|2|3]				#0=no binning, 1= do metaBat2 binning, 2=SemiBin, 3= MetaDecoder (experimental) (Default: 0)
	-BinnerMem [#]			define binning memory; e.g. 600 (Default: automatic)
	-BinnerCores [#]			cores used for Binning process (and checkM)
	-redoEmptyBins [0|1]		mostly for debugging: redo every sample where no bins where found (note that in some metagenomes there might be no bins)
	-checkM2 [0|1]			using checkM2 to assess bin quality
	-checkM1 [0|1]			using checkM1 to assess bin quality

# mapping
	-mapper [1|2|3|4]				1: bowtie2, 2:bwa, 3: minimap2, 4:kma, -1:auto (Default: -1)
	-mappingCores [#]				cores # used for mapping
	-mappingMem [#]				memory # used for mapping bwa/bwt2 in GB (Default: auto)
	-mapSortMem  [#]			memory # used for samtools sort in GB (Default: auto)
	-mappingCoverage			1: calculate coverage per predicted gene, contig, windows (Default: 1)
	-mapSupportReadsOntoAssembly [#]   1: map also support reads (e.g. PacBio in case of hybrid assemblies) onto final assembly
	-rmDuplicates [0|1]			1: remove read duplicates (Default: 1)
	-mapperFilterIll [# # #]				parameters for postprocessing mappings for short read data (Default: "0.05 0.75 20")
	-mapReadsOntoAssembly [0|1] 	1: map original reads to assembly to estimate contig (and bin, gene etc) abundance (Default: 1)
	-saveReadsNotMap2Assembly [0|1]			1: save reads not mapping to assembly in separate file (Default: 0)
	-remap2assembly [0|1]			1: redo the mapping to assembly (Default: 0)
	-JGIdepths [0|1]			1: calculate jgi coverage, only required when using MetaBAT2 binning (Default: automatic)

# SNPs
	-getAssemblConsSNP [0|1]		1: SNPs (onto self assembly); calculates consensus SNP of assembly (useful for checking assembly gets consensus and Assmbl_grps)
	-get2ndMappingConsSNP [0|1]		1: calculate consensus SNPs for mappings against references (map2tar mode)
	-redoAssmblConsSNP [0|1]		1: redo getAssemblConsSNP (Default: 0)
	-redoGeneExtrSNP [0\1]			1: redo gene extractions from consensus SNP contig (Default: 0)
	-SNPjobSsplit [#]				#: split consensus SNP jib further (Default: 1)
	-SNPsaveVCF [0|1]				1: save SNPs to VCF file (Default: 0)
	-SNPcaller [MPI|FB]				"MPI" mpileup or ".FB" for freebayes (Default: MPI)
	-SNPcores [#]				number of cores ‘#’ used for SNP calling
	-SNPmem [#]				Memory allocated for consensus SNP calling process in Gb (Default: 23)

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


### I/O (Input/Output): important considerations and design decisions

<details>
  <summary>Expand section</summary>

Analysing a shotgun metagenomic experiment can be a computationally extremely demanding task, as in some experiments several TB of data can be accumulating. MG-TK was designed with the latter case in mind, but can of course also handle smaller experiment.  
In order to be able to cope with these data amounts, a lot of 'file juggling' is happening behind the scenes. A lot of temporary files are being created that don't need to be saved on long term storage solution that are backed up and generally also slower. For this purpose big servers usually have a 'scratch' dir that is the global temporary storage on which all nodes in a cluster can write, but that is not backed up and might be cleaned infrequently. Further, usually each node has a local temp dir, to which only that specific node has access. Using these temporary solutions does make the whole cluster more stable and also enable other users to use a cluster more efficiently. To give you an example: if you have an IO heave process like searching with diamond through a lot of reads, you will use up the bandwith provided by your permanent storage very quickly. This could lead to situations where 500 cores on the cluster are busy with running in parallel diamond searches, but since the IO is so severely limited, only a small fraction of data trickles through to these jobs, effectively maybe giving 16 cores work. In this case the cluster would be unnecessarily blocked and the 500 core job would also take much longer than needed. That is the reason why file juggling is so important and why so much development effort went into optimizing this for MG-TK.  
To take advantage of this, I strongly recommend to ask your sysadmin where the local and global temp storage on your cluster are and set in the MG-TK config the variables 'globalTmpDir' and 'nodeTmpDir' variables correspondingly. 

</details>

### Useful configurations to track and check on MG-TK jobs

The most common reason why MG-TK jobs fail are related to node configurations (available ram, hdd space, CPUs). There are several alias' that are usful in checking on slurm jobs that are running on your local HPC, understanding how MG-TK processes your samples and fixing errors. Thus following up jobs and checking their error logs is essential in understanding limitations in your current environment and get your metagenomes processed effectively, as listed below:

<details>
  <summary>Expand section</summary>

These aliases can be directly added to your ~/.bashrc (just make sure the .bashrc is loaded):

```{sh}
#list running jobs with more relevant info
alias sq='squeue -u $USER -o "%8i %.4P %.14j %.2t %8M %.3C %.15R %20E"'
#check where job bash, std output, error output is stored, dependencies etc
alias si='scontrol show job'
#delete jobs that have DependencyNever status
alias scDN="squeue -u $USER | grep dencyNev | cut -f11 -d' ' | xargs  -t -i scancel {}"
#show the number of jobs currently running for different users on your cluster; useful for estimating how busy the HPC currently is
alias busy="squeue | sed -E 's/ +/\t/g' | cut -f5 | sort | uniq -c | sed -E 's/ +//' | sort -k1 -n -t' '"
#show output log of job
sio() {
JID=$1
if test "$#" -eq 0; then
JID=$(squeue -u hildebra | grep $USER | grep -v 'interact' | awk '{$1=$1};1' | cut -f1 -d' ' | head -1)
fi
cat $(scontrol show job $JID | grep 'StdOut' | sed 's/.*=//g')
}
#show error log of job
sie() {
JID=$1
if test "$#" -eq 0; then
JID=$(squeue -u hildebra | grep $USER | grep -v 'interact' | awk '{$1=$1};1' | cut -f1 -d' ' | head -1)
fi
cat $(scontrol show job $JID | grep 'StdErr' | sed 's/.*=//g')
}
#show bash script (commands) of job
sis() {
JID=$1
if test "$#" -eq 0; then
JID=$(squeue -u hildebra | grep $USER | grep -v 'interact' | awk '{$1=$1};1' | cut -f1 -d' ' | head -1)
fi
cat $(scontrol show job $JID | grep 'Command' | sed 's/.*=//g')
}
```
</details>

## Additional usage scenarios

MG-TK can be used for a bulk of tasks not directly related to initial assembly, profiling or MAGs, but often related to postprocessing these. Two usage scenarios (map2tar and building phylogenies) are listed below.

<details>
  <summary>Expand section</summary>

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

my $bts = "/path/to/MATAF3/secScripts/phylo/buildTree5.pl";
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

</details>


## Trouble shooting

<details>
  <summary>Expand trouble shooting section</summary>

This section lists a number of typically occurring problems that are usually not addressable by programming/bug fixing. Please look here first if an error you encountered is already listed.


### Setting environmental variables

**Problem:** If you have issue to source or define variables from/inside of MFF micromamba environment.

**Solution:** Adding 2 files into your MFF micromamba environment will help you resolve the issue (talk to Klara or Jogi).

### sbatch error: Batch job submission failed: Job dependency problem

**Problem:** Usually some jobs run but some jobs are hanging out in queue and never finish. In that case we can use `scontrol show job [ID]` to see which are the jobs and dependencies. 

**Solution:** If the dependencies are never fulfilled then we can delete all those jobs using `scancel`, after that rerun MG-TK. What MG-TK does, is to pick up where it left of - the pipeline is designed in a way that it creates `stone` files for finished processes so the pipeline knows where to continue. Sometimes files are also still in scratch dir and not copied to final dir. You just need to re-run command.

You can also delete **all** jobs where the dependency failed, saving you a lot of time (command: `squeue -u $USER | grep 'dencyNev' | cut -f11 -d' ' | xargs  -t -i scancel {}`) or ask MG-TK to do this in case the max number of jobs is reached via the flag `-killDepNever 1`.

### local tempdir on node not writable

**Problem:** Sometimes you will encounter an error where a job starts for ~1 sec on a node and immediately aborts, the error log (.etxt) showing something about not being able to create/write to a local dir (e.g. `/nbi/local/tmp/12312421/MG-TK`). This is usually the local SSD space not being available (for various reasons not related to MG-TK). Such nodes act as kind of 'honeypots', accepting a lot of jobs and killing them immediately. This can seriously harm your performance in getting jobs done.

**Solution:**  Therefore it is important to a) let you local sysadming know that the SSD is no longer available on said node (node name is always printed as first line in the MG-TK .otxt logs for a job) and b) you can exclude this node from MG-TK using the `-excludeNodes [nodename1,nodename2]` MG-TK flag.

### Recursion error while running contig stats

**Problem:** While calculating contig stats for large samples, you may encounter `RecursionError: maximum recursion depth exceeded in comparison` error. 

**Solution:** Increase the recursion depth in `extract_gtdb_mg.py` by using `sys.setrecursionlimit()` function, e.g. include `sys.setrecursionlimit(1500)` at the beginning of the python script. 

### Unusually high quality values

**Problem:** Some samples show unusually high quality of over 80

**Solution:** Some sequencing systems omit detailed quality values over a certain threshold. In this case, the quality values given are among the highest possible (>80) and do not represent the actual quality of sequences. To solve this quality values have to be calculated manually.

### GeneCat stops without producing MGS and no error messages appear

**Problem:** GeneCat.pl stops without error messages but no MGS are produced. Close examination shows missing genes among some bins, but these genes are present in the assembly.

**Solution:** This error can occur when previous GeneCat runs stop unexpectedly or fail due to previous issues with the assembly. In this case, files could be created but not completed. The MGS's are then not able to finish due to the expectation of some genes being present but the previous run had stopped before these could be written. It is best to make sure all assemblies are completed and then start a fresh GeneCat run.

### Automatic installation with the installer script does not finish due to several dependency issues

**Problem:** Parts of the installer using micromamba do not finish due to dependency issues

**Solution:** This can be resolved stepwise. First, make sure the environment where the problem occurs is created with the right name in micromamba. Then install packages (or dependencies) that have issues manually with micromamba. Try conda_forge first, then bioconda for the -c parameter. Restart the installer and note down any further issues. When a problem occurs with a package that is already installed, it can help to remove it and then reinstall it manually. If manual installation does not resolve conflicts, remove version numbers from the yml file of problematic packages and start this process again. It is important to note this somehow, in order to troubleshoot later on if any incompatible version of a package was installed this way.

</details>

## License, citations etc

### Used software
plenty.. please refer to helpers/install/MG-TK.yml, helpers/install/GTDBTK.yml, helpers/install/checkm2.yml and reqPackages.R

Other important software used:
sdm, LCA: both part of https://lotus2.earlham.ac.uk/

### Citing MG-TK

**Please cite MG-TK with:**
- Assembly mode: Hildebrand, F. et al. Antibiotics-induced monodominance of a novel gut bacterial order. Gut 68, 1781–1790 (2019). 
- Strain mode: Hildebrand, F. et al. Dispersal strategies shape persistence and evolution of human gut bacteria. Cell Host & Microbe 29, 1167-1176.e9 (2021). 
- Assembly-independent mode: Bahram, M. et al. Metagenomic assessment of the global diversity and distribution of bacteria and fungi. Environmental Microbiology 23, 316–326 (2021).
- sdm, LCA: Özkurt, E. et al. Microbiome (2022).

Falk Hildebrand <Falk.Hildebrand@gmail.com>

### License

 Copyright (c) 2017-2023 Falk Hildebrand

 MG-TK is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 2 of the License, or
 (at your option) any later version.

 MG-TK is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the file LICENSE for more details.

 You should have received a copy of the GNU General Public License
 along with the source code.  If not, see <http://www.gnu.org/licenses/>.
