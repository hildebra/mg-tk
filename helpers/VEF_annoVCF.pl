#!/usr/bin/perl
use strict;
use warnings;

my $gff = "/g/bork5/hildebra/results/TEC2/SebSNP/TEC2_RAST.gff";
my $fna = "/g/bork5/hildebra/results/TEC2/SebSNP/MM3.TEC2.scaffs.fna";
my $gbk = "/g/bork5/hildebra/results/TEC2/SebSNP/TEC2_RAST.gbk";

#my $inVCF= "/g/bork5/hildebra/results/TEC2/SebSNP/TEC2snpSEB.txt";
my $inVCF= "/g/bork5/hildebra/results/TEC2/SebSNP/MM3.fb.vcf";
my $outVCF= "/g/bork5/hildebra/results/TEC2/SebSNP/TEC2snpSEB_VEF.vcf";

my $outDir = "/g/bork5/hildebra/results/TEC2/SebSNP/coovar/";
system "mkdir -p $outDir";
my $gffgz = $gff.".gz";

#VEF
my $VEFdir = "/g/bork3/home/hildebra/bin/ensembl-tools-release-81/scripts/variant_effect_predictor/";
my $VEFbin = "$VEFdir./variant_effect_predictor.pl";
my $CACHbin = "$VEFdir./gtf2vep.pl";
my $CAHCdir = "$VEFdir/custCache/";
#system "mkdir -p $CAHCdir";
#die "sort -k1,1 -k2,2n -k3,3n $gff | /g/bork3/home/hildebra/bin/samtools-1.2/tabix-0.2.6/./bgzip > $gffgz";
#die "/g/bork3/home/hildebra/bin/samtools-1.2/tabix-0.2.6/tabix -p gff $gffgz";
#build custom cache
#die "$CACHbin -i $gff -f $fna -d 81 -s TEC2 --dir $CAHCdir";
#die "$VEFbin --custom $gff,myFeatures,gff,overlap,0 -i $inVCF -o $outVCF --force_overwrite --offline --species TEC2 --everything --format vcf --vcf --dir $CAHCdir";


#coovar: problems with creating indices
# die "/g/bork3/home/hildebra/bin/coovar-0.07/./coovar.pl -e $gff -r $fna -v $inVCF -o $outDir\n";
 
 #snpeff
 my $SEFFbin = "/g/bork3/home/hildebra/bin/snpEff/snpEff.jar";
 my $SEFFcfg = "/g/bork3/home/hildebra/bin/snpEff/snpEff.config";
 my $SEFFdata = "/g/bork3/home/hildebra/bin/snpEff/data/";
 my $geno = "TEC2.1"; my $Gname = "TEC2";
 my $custDir = "$SEFFdata/$geno/";
 if (0){ #build DB anew..
	my $cfgAdd = "# $Gname, version 1, test
	$geno.genome : $Gname ";
	if (`cat $SEFFcfg` =~ /$geno\.genome : $Gname /){
		print "Entry for such a genome already exists in DB"; exit(3);
	} else {
		#system "echo \"$cfgAdd\" >> $SEFFcfg";
	}
	system "mkdir -p $custDir; cp $fna $custDir/sequences.fa; cp $gff $custDir/genes.gff; cp $gbk $custDir/genes.gbk";
	my $cmd = "java -jar $SEFFbin build -genbank -v $geno\n"; 
 }
 my $cmd = "java -Xmx4g -jar $SEFFbin $geno $inVCF > $outVCF";
 die $cmd."\n";
 
 
 
 
 