#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -o 3.ill.GeneCat.sh.otxt
#SBATCH -e 3.ill.GeneCat.sh.etxt
#SBATCH --mem=81920
#SBATCH --export=ALL
#SBATCH -p "ei-medium,nbi-medium,ei-long,qib-medium,qib-long"
#SBATCH -J GjnGeCat
echo $HOSTNAME;
set -eo pipefail
ulimit -c 0;



echo ""
echo "Gene Catalogue creation example for illumina data"
echo "Relies on running and completing first: \"1.runMGTK_illumina.mfc\"  "
echo "Extension of simple example (nonsensical data) from script \"1.runMGTK_illumina.mfc\""
echo "By default, these results will be saved in $MGTKDIR/examples/output/ dir for this test."
echo "    However, this is not recommended for larger runs, where it is better to keep the input data and results dir separate from the MG-TK dir."
echo ""

if [ ! -e $MGTKDIR/examples/output/1.testAG/S4qiaS2/assemblies/metag/ContigStats/Coverage.pergene.gz ]; then
	echo "Test file not found (output/1.testAG/S4qiaS2/assemblies/metag/ContigStats/Coverage.pergene.gz)"
	echo "Ensure that 1.runMGTK_illumina.mfc finished correctly"
	exit
fi 

if [ ! -e $MGTKDIR/examples/output/1.testAG/S4qiaS3/assemblies/metag/Binning/SB/S4qiaS3.cm2 ]; then
	echo "Test file not found (output/1.testAG/S4qiaS3/assemblies/metag/Binning/SB/S4qiaS3.cm2)"
	echo "Ensure that 1.runMGTK_illumina.mfc finished correctly"
	exit
fi 



echo "It is recommended that you submit this run (sbatch 3.ill.GeneCat.sh), alternatively you could wait for the various submissions to finish"
echo ""
echo ""

#creates gene catalog in the specified outdir with specified cores, attempting to reuse existing dirs (in case catalog creation failed):

MAP=$MGTKDIR//examples/maps/testAG.map
GCoutdir=$MGTKDIR//examples////output///3.testGCill/


perl $MGTKDIR/secScripts/geneCat.pl -map $MAP -GCd $GCoutdir -mem 50 -cores 12 -clusterID 95 -doStrains 1 -continue 1 -binSpeciesMG 2 -useCheckM2 1 -useCheckM1 0 -MGset GTDB

