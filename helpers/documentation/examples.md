

## Examples

  
### MG-TK example dataset

We have prepared an example dataset that can be run directly after installing MG-TK and configuring it (see above). This example will 1) download public short and long read metagenomes 2) assemble short reads and 3) assemble short+long reads (hybrid assembly).

Please go to the directory mg-tk/examples/

To download all required data, run first 
```{sh}
bash 0.getExmplData.sh
```

After this is finished (check in the newly created mg-tk/examples/data/ dir for ~1.3Gb of data), you can either run 1.runMGTK_illumina.mfc (short read metagenomics) or 2.runMGTK_hybrid.mfc (short+long reads). Note that these are non-seniscal examples, i.e. the short and long reads are from completely independent experiments, don't expect interpretable results, this is purely to check if the technical process can run to completion.

How do you know everything finished as it should? Wait until all submitted jobs have finished, run the 1. or 2. script again until it reports that nothing is left to do. (Note:kill eventual "DependencyNeverSatisified" jobs for 1-2 times, if persists there might be a problem with runnning certain programs, where you need to start checking error logs, see Q&A below). Once everything is finished, MG-TK's postprocessing should complete and you will be able to find metagStats.txt and metagStatsReport.html - browse these files to get an overview of the data. 

In the next section we will give examples on how to create your own MG-TK runs.
  
  
### MG-TK metagenomic assembly and gene catalog

The figure below shows example of steps involved in the assembly-dependent mode. White rectangles indicate inputs and outputs, grey boxes name each of the steps, and yellow boxes show names of the scripts that are generated and submitted in each step. Blue boxes indicate additional steps that are required for subsequent MGS analysis.

<img src="assembly-dependent.svg" style="width: 800px;"/>


#### 1. create mapping file

Typically you would use Excel to create the mapping file and copy-paste it later into a text file (will be by default tab-delimited). This text file, typically with the file ending **.map** can then be saved to the HPC.

#### 2. make script with first command `RUN.mfc`

Insert your MG-TK command, a bash and slurm header in RUN.mfc. 

Example:

```{sh}
#!/bin/bash
#SBATCH -J SUB_MF
#SBATCH -N 1 --cpus-per-task=1 --mem=10024 --export=ALL
#SBATCH -o [currentDir]/run_mgtk_mhit.mfc.otxt
#SBATCH -e [currentDir]/run_mgtk_mhit.mfc.etxt
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

- Run with `bash run_mgtk_mhit.sh` --> this will submit a lot of different jobs to the HPC queue

 console output:
```
        This is MG-TK 0.33
        Using qsubsystem: slurm
        Using qsubsystem: slurm
        /projects/data/results/mgtk_test1/LOGandSUB/qsub.log
        Reset range of samples to 40

        ======= Mouse11T0 - 0 - M11T0 =======
        1:2  1:0
        SUB:_UZ0        SUB:_SDM0       SUB:_cln0       Running Contig Stats on assembly

        ======= Mouse11T1 - 1 - M11T1 =======
        2:2  1:0
        SUB:_UZ1        SUB:_SDM1       Assembly stepSUB:_A1    SUB:_GP1        SUB:_cln1       Running Contig Stats on assembly
        SUB:_CS1
```

- more advanced usage: run `sbatch run_mgtk_mhit.sh`. This will submit the job to the cluster queue, and from there the .mfc job will submit more jobs. The output from MG-TK will be stored in `#SBATCH -o [currentDir]/run_mgtk_mhit.mfc.otxt` and `#SBATCH -e [currentDir]/run_mgtk_mhit.mfc.etxt` defined above.


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
perl /hpc-home/project/mg-tk/secScripts/geneCat.pl \
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



### Assembly-independent MG-TK mode


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


### Hybrid assemblies

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

#### 2. setup MG-TK run

Setup your MG-TK run like you would normally setup an assembly dependent metagenomic analysis (see examples above). However, these flags should be defined:
 - `-mapSupportReadsOntoAssembly 1`: setting this to `1` will lead to support reads (PB reads in this case) being mapped onto the assembly. A coverage profile is created that is then **separately** used from the illumina coverage profile in the binning step (which in our experience can significantly boost the recovery of MAGs).
 - `-mapper -1`: set the choice of mapper to the default (`-1`) which means that a mapper will be automatically selected. You could set this to e.g. **1** to do all mapping with bowtie2, or **3** to use minimap2 everywhere, but **-1** is the recommended choice as this will use bowtie2 for ill reads and minimap2 for PB reads.
 - `-assembleMG 5`: this flag is crucial and tells MF to conduct a hybrid assembly, using **megahit for illumina** and **metaMDBG for PB** reads.
 - `-inputReadLengthSuppl 8000`: not crucial, but good to have. Here we estimated that our support PacBio reads are on average 8000 bp long.
 
Your example MG-TK call could now look like:
 
```{sh}
MAP=/path/to/map/PB_hybrid.map

perl $MGTKDIR/MG-TK.pl -map $MAP -inputFQregexSingle '.*\.fastq\.gz' -inputFQregex1 '(.*_R1_001\.fastq\.gz)|(.*[_\.]1\.f[^\.]*q\.gz)$' -inputFQregex2 '(.*_R2_001\.fastq\.gz)|(.*[_\.]2\.f[^\.]*q\.gz)$' -inputBAMregex '.*\.bam$' \
-assembleMG 5 -spadesCores 12 -spadesKmers "25,43,67,87,111,127" -spadesMemory 200 -MetaBat2 2 \
-mapper -1 -mapSupportReadsOntoAssembly 1 \
-filterHostKrak2DB /path/to/kraken2/hsap/ -filterHostRds 1 \
-getAssemblConsSNP 1 -rmSmplLocks 0 \
-submit 1 -inputReadLength 150 -inputReadLengthSuppl 8000 \
-from 0 -to 1
```

#### 3. running MG-TK in hybrid mode

We opted for a bit of a complicated processing of hybrid assemblies, that in the end allows for both usage of existing paths in MG-TK as well as supporting more complicated (assemblyGrps) sample setup.
What this means for you as the user is simply that you need to run the above MG-TK command several times. In the first iteration this will trigger the megahit illumina-assembly, in the second iteration remaining samples in assemblyGrp are mapped onto this illumina-assembly, in the third iteration the hybrid assembly with metaMDBG is started, fourth iteration will map remaining samples in assemblyGrp onto hybrid-assembly, fifth iteration will then finally starting consensus SNP calling and binning. Remember to wait between iterations until all current jobs have finished (though sample locks should normally prevent double submissions).

_tldr; this mode requires several iterations to complete_ 


