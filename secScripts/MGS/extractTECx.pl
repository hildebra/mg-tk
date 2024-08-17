#!/usr/bin/perl
#script to extract all the TEC2 from a new gene catalog & describes subsequent work
#run TAMOC
#./TAMOC.pl
#run GC
#./geneCat.pl /g/bork5/hildebra/data/metaGgutEMBL/MM_at_v3_T2subset.txt,/g/bork5/hildebra/data/metaGgutEMBL/ABRtime.txt,/g/bork5/hildebra/data/metaGgutEMBL/T2_HMP.txt /g/scb/bork/hildebra/SNP/GCs/T2_HM3_GNM3_ABR 1 95
#run Canopy clustering
#./geneCat.pl MGS /g/scb/bork/hildebra/SNP/GCs/
#extract marker genes
use strict; use warnings;
use threads ('yield',
                 'stack_size' => 64*4096,
                 'exit' => 'threads_only',
                 'stringify');




#binary configuration
my $MATAdir = "/g/bork3/home/hildebra/dev/Perl/reAssemble2Spec/";
my $hdir = "$MATAdir/helpers/";
my $GenomeDir = "/g/bork3/home/hildebra/results/prelimGenomes/v4/";
my $onlyRfilt = 0;

system "mkdir -p $GenomeDir" unless (-d $GenomeDir);


if (0){
	#select MaxBin clusters, based on mOTU genes
	#"$hdir/GC/getMOTU_MG.pl Cluster1088";
	#/g/bork3/home/hildebra/results/TEC2/Shini_motu_TS/MGs//Cluster1088.MG.fna
}

#configuration of paths & samples
#this assumes that a set of genomes from MaxBin was already selected
#2=motu_linkage_group_838, 6=Cluster1088, 7=Cluster1104, 8=Cluster1630, 9=Cluster1576, 10=Cluster1092
my @inDs = ("alien-11-376-0/","alien-11-374-0/","alien-11-374-0/","alien-11-374-0/","alien-11-380-0/","alien-9-7-0/","alien-11-60-0","alien-11-898-0/","alien-11-877-0/");
my @SmplN = ("MM4","MM3","MM3","","MM14","MM35","MM11","MM18");#exclusively determines from which dir the scaffolding is made..
my @motuName = ("motu_linkage_group_838","","","","Cluster1088","Cluster1104","Cluster1630","Cluster1576","Cluster1092");
my @num=("2","3","4","5","6","7","8","9","10");
my @MBnum = ("2","2","3","4","2","2","2","2","3");
my $smplDir = "/g/scb/bork/hildebra/SNP/GNMass3/";
my $odir = "/g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/";
my $gcat= "/g/scb/bork/hildebra/SNP/GCs/T2_HM3_GNM3_ABR";
my $maps = "/g/bork5/hildebra/data/metaGgutEMBL/MM_at_v5_T2subset.txt";
my $LpreClus= "$gcat/LuisSpecs/";
my @thrs;#array to save threads in..
#helpers/./getMarkersMultiSmpl.pl FMGcrossSmpls alien-11-376-0/ 2 TEC2 /g/scb/bork/hildebra/SNP/GNMass3/ /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/ /g/scb/bork/hildebra/SNP/GCs/T2_HM3_GNM3_ABR 90
#helpers/./MGSonGeneFetch.pl /g/scb/bork/hildebra/SNP/GCs/T2_HM3_GNM3_ABR /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v4/T2/TEC2Names.txt

#get initial set of genes, based on maxbin and canopy clustering results
my $st = 4; my $en = @inDs;
for (my $j=$st;$j<$en;$j++){
	last if ($onlyRfilt);
	my $i = $num[$j];
	system "mkdir -p $odir/T$i/";
	my $cmd = "set -e\n";
	
	#get marker genes and their abundances from maxbin pre-genome (this can also be used to build a tree from assembly)
	#note -> could also be used to get per sample gene consensus seq for all genes associated to a species.
	$cmd .= "$hdir/./getMarkersMultiSmpl.pl FMGcrossSmpls $inDs[$j] $MBnum[$j] TEC$i $smplDir $odir/T$i/ $gcat 90 $maps > $odir/T$i/gMMS.log\n";
	#die $cmd ."\n";
	#next, get the gene matrix of maxbin and canopy clustering pre-genomes, and create gene abundance matrix for these
	#file $odir/T$i/TEC${i}Names.txt contains 40 (or more) FMG IDs
	#motuclus refers to binnings from luis, based on gene to mOTU abundance
	my $MOTUclus = "X"; $MOTUclus = "$LpreClus/$motuName[$j].tsv" if ($motuName[$j] ne "");
	print "T$i\t$MOTUclus\n";
	$cmd .= "$hdir/./MGSonGeneFetch.pl $gcat $odir/T$i/TEC${i}Names.txt  > $odir/T$i/MGS_GF.log\n";
	#
	$cmd .= "$hdir/./motu_genome_integrate.pl $gcat $motuName[$j] $MOTUclus > $odir/T$i/motu_genome.log\n";
	#die $cmd."\n";
	$thrs[$j] = threads->create(sub{system $cmd;});
	#print "$j ";
}
print "\n";
for (my $t=$st;$t<$en;$t++){
	last if ($onlyRfilt);
	my $state = $thrs[$t]->join();
	if ($state){print "Thread $t exited with state $state\n";}
}
exit();

#------------------------------------------------------
# after this run R filter (copy matrix files) and insert them in odir/TX/R_filt
#Final Genome assembly step: a) discovery or b) scaffolding
#------------------------------------------------------
my $exploratory = 1; my $doScaffolding = 0; my $doMinimus = 0;
for (my $j=$st;$j<$en;$j++){
	#last;
	my $i = $num[$j];
	#exploratory route
	my $cmd = "";
	if ($exploratory){
		#this part takes the R filtered gene list and extracts the corresponding contigs (for each sample)
		system "mkdir -p $odir/T$i/R_filt" unless (-d "$odir/T$i/R_filt");
		$cmd = "$hdir./geneListSameAssembly.pl stat $odir/T$i/R_filt/T${i}_filt.txt $gcat $odir/T$i/R_filt/contigs2 > $odir/T$i/R_filt/extractRun.txt\n";# $SmplN[$j]\n";
	} elsif ($doScaffolding) { #route for scaff based assemblies, where we know already which sample had the best scaff statistics
		#$cmd = "$hdir./geneListSameAssembly.pl stat $odir/T$i/R_filt/T${i}_filt.txt $gcat $odir/T$i/R_filt/contigs2 $SmplN[$j] > $odir/T$i/R_filt/extractSpecificsRun.log\n";# \n";
		#$cmd= "./TAMOC.pl scaffold /g/scb/bork/hildebra/SNP/GNMass3/TECtime/v5/T2/R_filt/contigs/MM4.ctgs.fna TEC$i MM4";
		next if ($SmplN[$j] eq "");
		#scaffold contigs (requires mate pair library)
		$cmd= "$MATAdir/./TAMOC.pl scaffold $odir/T$i/R_filt/contigs2/$SmplN[$j].ctgs.fna TEC$i $SmplN[$j]";
	} elsif ($doMinimus){ #minimus assembly across different samples: doesn't work, don't use
		#currently needs to be handedited
		$cmd= "$MATAdir/helper/miniSmpl.pl";
		die $cmd;
	}
	
	print $cmd."\n";
	system $cmd;
	#die;
}

die;

#prep scaffolds to map against
#just copy them to new dir, not complicated..

#my @chFs = ("MM3","MM3","","MM3","MM7");#v4
my @chFs = ("MM4","","","","");#v5
my @scArr ;my @TecN;
for (my $j=$st;$j<$en;$j++){ 
	my $i = $num[$j];
	system "rm -r $GenomeDir/T$i";
	system "mkdir -p $GenomeDir/T$i/" unless (-d "$GenomeDir/T$i/");
	my $tarFil = "$odir/T$i/R_filt/contigs/$chFs[$j].ctgs.fna";
	$tarFil = "/g/scb/bork/hildebra/SNP/GNMass3/alien-11-377-0/scaffolds/TEC$i/BESST_output/pass2/Scaffolds_pass2.fa";
	die "can't find tar file $tarFil\n" unless (-e $tarFil);
	my $genRestF = "$GenomeDir/T$i/TEC$i.ctgs.fna";
	system "cp $tarFil $genRestF";
	push(@scArr,"$genRestF");
	push @TecN,"TEC".$i;
}



#------------------------------------------------------
#Map reads to genome (decoy mapping)
#required for SNP calling process
# TODO: update command
#------------------------------------------------------

#./TAMOC.pl map2tar /g/bork3/home/hildebra/results/prelimGenomes/TEC3/MM3.TEC3.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC3/TEC3ref.fasta,/g/bork3/home/hildebra/results/prelimGenomes/TEC4/MM3.TEC4.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC4/TEC4ref.fasta,/g/bork3/home/hildebra/results/prelimGenomes/TEC5/MM4.TEC5.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC5/TEC5ref.fasta,/g/bork3/home/hildebra/results/prelimGenomes/TEC6/MM29.TEC6.scaffs.fna,/g/bork3/home/hildebra/results/prelimGenomes/TEC6/TEC6ref.fasta TEC3,TEC3r,TEC4,TEC4r,TEC5,TEC5r,TEC6,TEC6r
#/g/bork5/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.fa
print "./TAMOC.pl map2tar ".join(",",@scArr) ." ". join(",",@TecN)."\n";

#./TAMOC.pl map2tar /g/bork3/home/hildebra/results/TEC2/v5/TEC2.MM4.BEE.GF.rn.fa,/g/bork5/hildebra/results/TEC2/v5/T3/T3.mini2.3smpl.fna,/g/bork5/hildebra/results/TEC2/v5/T4/T4.mini2.3smpl.fna,/g/bork3/home/hildebra/results/TEC2/v5/T6/TEC6.ctgs.rn.fna T2d,T3d,T4d,T6d


#afterwards follows SNP calls
#./SNPcalls.pl
#and tree building
#helpers/./buildTree.pl
#helpers/./growthRate.pl




