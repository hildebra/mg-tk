#!/usr/bin/bash
#download example data for MG-TK
#(c) Falk Hildebrand

if [ -z ${MGTKDIR+x} ]; then echo "MGTKDIR is unset"; else echo "MGTKDIR is set to '$MGTKDIR', downloading data to $MGTKDIR/examples/data/"; fi
rm -rf $MGTKDIR/examples/data
mkdir $MGTKDIR/examples/data








echo ""
echo ""
echo "downloading PacBio example metagenomes"
echo ""

#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/009/SRR15489009/SRR15489009_subreads.fastq.gz --output $MGTKDIR/examples/data/SRR15489009_sub.fastq.gz
#zcat  $MGTKDIR/examples/data/SRR15489009_sub.fastq.gz | head -n 100000 >  $MGTKDIR/examples/data/SRR15489009_sub.fq.gz


curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/009/SRR15489009/SRR15489009_subreads.fastq.gz | zcat | head -n 300000 | gzip > $MGTKDIR/examples/data/SRR15489009_sub.fq.gz

#curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/013/SRR15489013/SRR15489013_subreads.fastq.gz --output $MGTKDIR/examples/data/SRR15489013_sub.fastq.gz
#zcat  $MGTKDIR/examples/data/SRR15489013_sub.fastq.gz | head -n 100000 >  $MGTKDIR/examples/data/SRR15489013_sub.fq.gz

curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR154/013/SRR15489013/SRR15489013_subreads.fastq.gz | zcat | head -n 300000 | gzip > $MGTKDIR/examples/data/SRR15489013_sub.fq.gz


#rm -f $MGTKDIR/examples/data/SRR15489013_sub.fastq.gz $MGTKDIR/examples/data/SRR15489009_sub.fastq.gz




echo ""
echo "downloading illumina example metagenomes (PRJNA529586)"
echo ""
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/002/SRR8797712/SRR8797712_1.fastq.gz --output $MGTKDIR/examples/data/SRR8797712_1.fq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/002/SRR8797712/SRR8797712_2.fastq.gz --output $MGTKDIR/examples/data/SRR8797712_2.fq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/003/SRR8797713/SRR8797713_1.fastq.gz --output $MGTKDIR/examples/data/SRR8797713_1.fq.gz
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR879/003/SRR8797713/SRR8797713_2.fastq.gz --output $MGTKDIR/examples/data/SRR8797713_2.fq.gz
 
 
