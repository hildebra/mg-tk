#!/usr/bin/env perl
#annotates specI's in the dataset
#also creates specI abundance tables
#will create tax abundance table 
#dont forget to add custom genomes to specI database!! -> this was now replaced by 3rd arg being a MGS list
#  perl annotateMGwSpecIs.pl /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR 12 /g/bork3/home/hildebra/data/SNP/GCs/T2_HM3_GNM3_ABR/Canopy3/clusters.txt
# perl annotateMGwSpecIs.pl /g/bork3/home/hildebra/data/SNP/GCs/alienGC2 12 /g/bork3/home/hildebra/data/SNP/GCs/alienGC2/Binning/MetaBat/MB2.clusters.core /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/extended_Tax.txt
# perl annotateMGwSpecIs.pl /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/ 12 /g/bork3/home/hildebra/data/SNP/GCs/DramaGCv5/Binning/MetaBat/MB2.clusters.ext.can.Rhcl.mgs.srt /g/scb/bork/hildebra/SNP/GCs/DramaGCv5/Binning/MetaBat/extended_Tax.txt


use warnings;
use strict;

use Mods::GenoMetaAss qw( readClstrRev systemW median readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt );
use Mods::geneCat qw(calculate_spearman_correlation read_matrix correlation);

use Mods::IO_Tamoc_progs qw(getProgPaths);
use List::Util;


sub readMotuTax; sub fixGTDBtax;
sub readGene2mlinkage;
sub lambdaBl;
sub readNCBI;sub read_speci_tax;
sub createProfileMGS;
sub MGSassign; sub transferSI2MGS;

#sub calculate_spearman_correlation;
#sub read_matrix; 

sub getCorrs;
sub passBlast;
sub add2geneList; sub rm4geneList;
sub sanityCheckCorr;
sub specImatrix;
sub rmFromList;
sub createAreadSpecItax; sub writeSpecItax;
sub disentangleMultiAssigns;

#my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI/"; my $freeze11=1;
my $freeze11=0;
my $doGTDBtax = 1; #Firmicutes_A etc
#progenomes.specIv2_2
my $globalCorrThreshold = 0.8; # determines cutoff, when still to accept correlating genes into specI
my $reblast=0;#do blast again?

my $SpecID=getProgPaths("specIPath0");#"/g/bork3/home/hildebra/DB/MarkerG/specI_f11/";
my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";
my $lambdaBin = getProgPaths("lambda");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda";
my $lambdaIdxBin = $lambdaBin."_indexer";#getProgPaths("");#"/g/bork3/home/hildebra/dev/lotus//bin//lambda/lambda_indexer";
my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $bts = getProgPaths("buildTree_scr");

if (@ARGV == 0){
	die "Not enough input args: use ./annotateMGwMotus.pl [path to GC] [# Cores]\n";
}

my $GCd = $ARGV[0]."/";
my $BlastCores = $ARGV[1];
my $MGSfile = ""; #MGS annotations; will take precedence of blast specI annotations..
$MGSfile = $ARGV[2] if (@ARGV > 2);
my $MGStax = ""; #MGS taxonomy (fitting to MGSfile).. only if given will try to calc abundance table
$MGStax = $ARGV[3] if (@ARGV > 3);
my $MGdir = "$GCd/FMG/";
my $outD = "$GCd/Anno/Tax/SpecI";
if ($MGSfile ne ""){$outD = "$GCd/Anno/Tax/SpecI_MGS/";}

system "mkdir -p $MGdir/tax" unless (-d "$MGdir/tax");
system "mkdir -p $outD" unless (-d "$outD");
#die $outD;
my $motuDir = "";#"/g/bork3/home/hildebra/DB/MarkerG/mOTU";
#load motu DBs...
#my ($hr1,$hr2) = readMotuTax("$motuDir/mOTU.v1.1.padded.motu.linkage.map");
#my %LG2motu = %{$hr1}; my %motu2tax = %{$hr2};
#$hr1 = readGene2mlinkage("$motuDir/mOTU.v1.1.padded.motu.map");
#my %gene2LG = %{$hr1};
#$hr1 = readNCBI("/g/bork3/home/hildebra/DB/NCBI/ncbi_tax_table_synonyms_2014-01-23.txt");
#my %NTax = %{$hr1};

#$hr1 = read_speci_tax("$SpecID/specI.tax");
#my %specItax = %{$hr1};
#annotate against DB using lambda

if (0){ #too general
	my $tar = "$GCd/compl.incompl.95.fna"; my $DB = "$motuDir/263MetaRef10MGv9.cal.v2.nr.padded.fna";	my $taxblastf = "$GCd/compl.incompl.95.motuAss.tmp.m8";
	lambdaBl($tar,$DB,$taxblastf);
}

my %specIid;my %specItax;
#open I,"<$SpecID/progenomes.specIv2_2";

#open I,"<$SpecID/prok-refdb-v11.0.0_specI-v2_clusters-v1.map" or die "Can't open in specI map\n";
my $inSImap = "$SpecID/freeze12_2_representatives_augmentedWhitelist_mapping.txt";
open I,"<$inSImap" or die "Can't open in specI map $inSImap\n";
while (<I>){next if (m/^#/);chomp; my @xx = split /\t/;$xx[1] =~s/,//g; $specIid{$xx[1]} = $xx[0];$xx[1]=~m/^(\d+)\./; $specItax{$xx[0]}=$1;}
close I;
my %specItaxM; #real matrix with tax levels
#die "$specItax{specI_v2_Cluster1309}\n";
my $specIfullTaxHR = createAreadSpecItax(\%specItax,"$SpecID/specI.tax3","$SpecID/bacTaxGTDB.tsv");
#writeSpecItax($specIfullTaxHR,"$SpecID/specI.tax3");

my %specIfullTax = %{$specIfullTaxHR};
my $xtrLab= "";$xtrLab= ".rep" if ($freeze11);

#assign each COG separately
system "mkdir -p $MGdir" unless (-d $MGdir);
my %FMGcutoffs = (COG0012=>94.8,COG0016=>95.8,COG0018=>94.2,COG0172=>94.4,COG0215=>95.4,COG0495=>96.4,COG0525=>95.3,COG0533=>93.1,COG0541=>96.1,
COG0552=>94.5,COG0048=>98.4,COG0049=>98.7,COG0052=>97.2,COG0080=>98.6,COG0081=>98,COG0085=>97,COG0087=>99,COG0088=>99,COG0090=>98.8,COG0091=>99,
COG0092=>99,COG0093=>99,COG0094=>99,COG0096=>98.6,COG0097=>98.4,COG0098=>98.7,COG0099=>98.9,COG0100=>99,COG0102=>99,COG0103=>98.4,
COG0124=>94.5,COG0184=>98.2,COG0185=>99,COG0186=>99,COG0197=>99,COG0200=>98.4,COG0201=>97.2,COG0202=>98.4,COG0256=>99,COG0522=>98.6);
my %specItaxname; my %SpecIgenes; 
my %COGDBLass; my %COGass; my %gene2COG;
my %gen2SIscore;

my %Q2S; #assignment of GC genes to specI id (poss. several ids)

my %COG2FMG; my %FMGlist;
open IC,"<$GCd/FMG.subset.cats" or die "Can't open $GCd/FMG.subset.cats\n";
while (<IC>){
	chomp;	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my @genes = split(/,/,$spl[2]);
	my $curCOG = $spl[0];
	$COG2FMG{$curCOG} = \@genes;
	foreach (@genes){$FMGlist{$_} = $curCOG;}
}
close IC;

#MGS related containers
my %Gene2MGS; my %speci2MGS; my %MGSlist;
my %speci2MGScnt; #counts how many different genes were used in total
if ($MGSfile ne ""){
	my %tmpTax;
	unless ($MGStax eq ""){
		open IT,"<$MGStax" or die "can't open MGS tax file $MGStax\n";
		while (<IT>){
			chomp;my @spl=split/\t/;
			next if ($spl[0] eq "domain" || $spl[0] eq "user_genome");
			my $id = shift @spl;
			if (@spl == 1){
				my $tmp=$spl[0]; $tmp =~ s/;;/;\?;/;$tmp =~ s/;$/;\?/;
				@spl = split /;/,$tmp;
				push(@spl,"?") while (@spl < 7);
			}
			@spl = fixGTDBtax(@spl) if ($doGTDBtax);
			$tmpTax{$id} = \@spl;
		}
		close IT;
	}
	my $taxF=0; my $taxN=0;
	open IM,"<$MGSfile" or die "Can't open MGS $MGSfile\n";
	while (<IM>){ 
		chomp; my @spl = split /\t/; my @spl2= split /,/,$spl[1]; 
		my $curMGS = $spl[0];
		foreach my $gen (@spl2){
			#chomp $gen;
			next unless (exists($FMGlist{$gen}));
			#print "$gen " if ($curMGS eq "MB2bin186");
			$Gene2MGS{$gen} = $curMGS;
			$MGSlist{$curMGS} = 1;
			#just insert in specI related objects.. just pretend it's great
			#die "$gen $curMGS\n" if ($gen <5);
			push(@{$SpecIgenes{$curMGS}{$FMGlist{$gen}}},$gen);
			$gen2SIscore{$gen}{$curMGS}=200;
		}
		#just gie it empty tax for now..
		if (!exists($specIfullTax{$curMGS})){
			if (exists($tmpTax{$curMGS})){
				$specIfullTax{$curMGS} = $tmpTax{$curMGS};
				$taxF++;
			} else {
				$specIfullTax{$curMGS} = ["Bins","?","?","?","?","?","?"];
				$taxN++;
			}
		}
	}
	close IM;
	print "Found ". scalar(keys %MGSlist)." MGS with $taxF/".($taxF+$taxN) ." taxonomies.\n";
}
undef %FMGlist;
#die;
#my @MGS = keys %Gene2MGS;die @MGS."\n$Gene2MGS{$MGS[0]}\n";
my $allOK=1;
foreach my $COG (keys %COG2FMG){#(@catsPre){
	my %specIcnt;
	my @genes = @{$COG2FMG{$COG}};
	#actual blast (heavy)
	my $cmd = "$samBin faidx $GCd/compl.incompl.95.fna ". join (" ", @genes ) . " > $MGdir/$COG.fna\n";
	system $cmd unless (-e "$MGdir/$COG.fna");
	$cmd = "$samBin faidx $GCd/compl.incompl.95.prot.faa ". join (" ", @genes ) . " > $MGdir/$COG.faa";
	system $cmd unless (-e "$MGdir/$COG.faa");	
	
	#print "$MGdir/tax/$COG${xtrLab}.tmp.m8\n";
	lambdaBl("$MGdir/$COG.fna","$SpecID/$COG${xtrLab}.fna","$MGdir/tax/$COG${xtrLab}.tmp.m8"); #.rep
	#next;
	
	die "no cutoff for $COG\n" unless (exists $FMGcutoffs{$COG});
	my $reqID = $FMGcutoffs{$COG};
	my $prevID="";my %finiM=();my $finiN=0;my $maxSI="";
	my $inM8 = "$MGdir/tax/$COG${xtrLab}.tmp.m8"; my $problemFlag=0;
	open I,"<$inM8" or die $!;
	while (my $line = <I>){
		chomp $line; 
		my @spl = split /\t/,$line;
		my $gen = $spl[0];
		if ($prevID ne $gen){ #reset some stats counters for MGS
			$prevID = $gen; %finiM=();$maxSI=""; $finiN=0;
		}
		$spl[1] =~ m/^([^\.]*\..[^\.]*)\./;
		die "can't find specI $1\n$inM8\n" unless (exists $specIid{$1});
		#get specI assignments
		my $speci= $specIid{$1};
		#only count each speci per gene once.. aim is to find possible candidates
		#TODO: could lead to too many hits and confusion..
		next if (exists($finiM{$speci}));# && $maxSI ne $speci);
		#next unless ($spl[1] =~ m/6666666/);
		my $pBval = passBlast(\@spl,$reqID,$inM8);
		if ($pBval == -2){
			close I;$allOK=0; print "Problem!\n";
			system "rm $inM8";$problemFlag=1;
			#rerun lambda..
			lambdaBl("$MGdir/$COG.fna","$SpecID/$COG${xtrLab}.fna","$MGdir/tax/$COG${xtrLab}.tmp.m8"); 
			last;
		}
		next if ($pBval==0);
		#next if (exists $Q2S{$spl[0]});
		$gene2COG{$gen} = $COG;
		$gen2SIscore{$gen}{$speci} = $spl[2];#perID

		#just link the specI to MGS..
		if (exists($Gene2MGS{$gen})){
			my $curMGS = $Gene2MGS{$gen};
			$speci2MGS{$curMGS}{$speci}++;
			if ($finiN == 0){$speci2MGScnt{$curMGS} ++ ; $finiN=1;}
			$finiM{$speci}=1;
			#next;
		}
		
		#save all valid hits  (%id) to a given specI
		if (!exists($Q2S{$gen}{$speci} )){
			$Q2S{$gen}{$speci} = scalar(keys(%{$Q2S{$gen}})); #print "m";
		}
		
		if (exists($SpecIgenes{$speci}{$COG} )){
			$COGDBLass{$COG}++;
		} 
		push(@{$SpecIgenes{$speci}{$COG}},$gen);
	}
	close I;
}

if ($allOK==0){die"One or more Blast files were not ok.. restart procedure\n";}

#write out SpecI -> MGS assignments (and tax)
my $SI2MGShr = MGSassign();
#transfer gene assignments to MGS, transfer SI tax -> MGS tax
transferSI2MGS($SI2MGShr);
my %SI2MGS = %{$SI2MGShr};


#die;
#print "Done initial blast\n$MGdir\n";
#foreach (sort {$COGass{$a} cmp $COGass{$b}} keys %COGass){print "$_ $FMGcutoffs{$_}:  $COGass{$_}($COGDBLass{$_})\n";}
my $hr1 = read_matrix("$GCd/Matrix.FMG.mat"); #here I need the coverage matrix..
#my $hr1 = read_matrix("$GCd/Mat.cov.FMG.mat"); #here I need the coverage matrix..
#my $hr1 = read_matrix("$GCd/Mat.med.FMG.mat"); #here I need the coverage matrix..
my %FMGmatrix= %{$hr1};

my %gene2specI; my %specIprofiles; my %SpecIgenes2;
#sort out best multi hit by correlation analysis
getCorrs(\%FMGmatrix,\%SpecIgenes);

disentangleMultiAssigns();
#my $T2cnt=0;foreach my $k (keys %gene2specI){$T2cnt++ if ($gene2specI{$k} eq "TEC2");}print "T2cnt = $T2cnt\n";

#check that all corrs here check out well..
sanityCheckCorr();
#my $T2cnt=0;foreach my $k (keys %gene2specI){$T2cnt++ if ($gene2specI{$k} eq "TEC2");}print "T2cnt = $T2cnt\n";

#second round.. go through blast again and reassign everything, that hasn't got a hit yet

undef %Q2S; my $xtraEntry=0;
#@catsPre = split/\n/,`cat $GCd/FMG.subset.cats`;
foreach my $COG (keys %COG2FMG){#(@catsPre){
	my %specIcnt;
	my @genes = @{$COG2FMG{$COG}};
	my $reqID = $FMGcutoffs{$COG};
	#secondary scanning of assignments...
	my $inM8 = "$MGdir/tax/$COG${xtrLab}.tmp.m8";
	open I,"<$inM8" or die "Can't open $inM8\n";
	while (my $line = <I>){
		chomp $line; 
		my @spl = split /\t/,$line;
		my $gid = $spl[0];
		#1 gene not already assigned in first high confidence pass
		next if (exists $gene2specI{$gid});
		#2 blast good enough
		next if (passBlast(\@spl,$reqID,$inM8)==0);
		
		$spl[1] =~ m/^([^\.]*\.[^\.]*)\./;
		die "can't find specI $1\n" unless (exists $specIid{$1});
		#get specI assignments
		my $speci= $specIid{$1};
		#direct translation to MGS, if SI correponds to this
		$speci = $SI2MGS{$speci} if (exists($SI2MGS{$speci}));
		
		
		#3 this MG has already been assigned in the high confidence initial assignments
		if (exists($SpecIgenes2{$speci}{$COG} )){next;}
		
		#4 correlate to species core, to make sure the gene kind of makes sense..
		if (exists($specIprofiles{$speci})){
			my $corr = correlation($specIprofiles{$speci},$FMGmatrix{ $gid }) ;
			next if ($corr < $globalCorrThreshold);
			#print $corr."\t";
		}

		if (!exists($Q2S{$gid}{$speci} )){
			$Q2S{$gid}{$speci} = scalar(keys(%{$Q2S{$gid}})); #print "m";
		}

		$specIcnt{$speci}++;
		$xtraEntry++;
		#and block gene slot
		add2geneList($speci,$COG,$gid);
		#$gene2specI{$spl[0]} = $speci;
	}
	close I;
	#print "Found ".keys(%Q2S)." assignments (".@genes.")\n";
	
	#deactivated..
	foreach my $gen (keys %Q2S){
		last;
		my @spl = sort {$Q2S{$gen}{$a} <=> $Q2S{$gen}{$b}} keys(%{$Q2S{$gen}});
		my $finalSpeci = $spl[0];
		my $set=0;
		#routine to distribute genes as good as possible between specIs
		for (my $i=(@spl-1);$i>=0;$i--){
			my $s = $spl[$i];
			if ( ($specIcnt{$s}>1 && $i!=0 ) || $set){#delete entries where too many specIs are available
				$specIcnt{$s}--;
			} else {
				$specIcnt{$s}--;
				$finalSpeci = $spl[$i]; $set=1;
			}
		}
		
		$gene2specI{$gen} = {$finalSpeci => undef};
	}
}

#die;
sanityCheckCorr();
#my @T2g = keys %{$SpecIgenes{TEC2}};print "T2: @T2g\n";
#my $T2cnt=0;foreach my $k (keys %gene2specI){$T2cnt++ if ($gene2specI{$k} eq "TEC2");}print "T2cnt = $T2cnt\n";


#write specI assignments for markerG
open O,">$MGdir/gene2specI.txt";
foreach my $k (keys %gene2specI){
	print O "$k\t". join(",",keys %{$gene2specI{$k}}) . "\t$gene2COG{$k}\n";
}
close O;
#add MGS profiles extra to the 
createProfileMGS();
#print @{$specIprofiles{MB2bin61}}."\n";
#create abundance profile
specImatrix("$outD/specI.mat",\%specIfullTax);

my %specIcnts; #just use for histgram
foreach my $k (keys %gene2specI){$specIcnts{ (keys %{$gene2specI{$k}})[0] }++;}
my %histo;
for my $k (sort {$specIcnts{$a} <=> $specIcnts{$b}} keys %specIcnts) {
	$histo{int $specIcnts{$k}/10}++;
   # print "$k $specIcnts{$k} $NTax{$specItax{$k}}\n" ;#if ($specIcnts{$k}>=40);   # bbb c aaaa
}
foreach (sort (keys %histo)){
	print "$_\t$histo{$_}\n";
}
print "Extra: $xtraEntry\n";

print "Finished SpecI annotations & matrix\n$outD\n";
exit(0);




#####################################################

sub MGSassign{

	my %Si2MGS ;
	return \%Si2MGS if ($MGSfile eq "");
	open OM,">$outD/MGS2speci.txt";
	
	foreach my $MGS (sort {$speci2MGScnt{$b} <=> $speci2MGScnt{$a}} keys %speci2MGScnt){
		print OM "$MGS";
		my @SIsAdd;
		my $valSI="";my $valStr=0;
		foreach my $si (sort {$speci2MGS{$MGS}{$b} <=> $speci2MGS{$MGS}{$a}} keys %{$speci2MGS{$MGS}}){
			 push(@SIsAdd,"$si:$speci2MGS{$MGS}{$si}");
			 if ($valSI eq "" && $speci2MGS{$MGS}{$si}>5){
				$valSI = $si;$valStr =$speci2MGS{$MGS}{$si}; 
			 } elsif ($valStr < 0.5*$speci2MGS{$MGS}{$si}){ #just too ambigous..
				$valSI = "";
			 }
		}
		if ($valSI ne ""){
			#print for file the link
			print OM "\t".join(";",@{$specIfullTax{$valSI}});
		} else {
			print OM "\t?;?;?;?;?;?;?";
		}
		print OM "\t".join(",",@SIsAdd)."\n";
		$Si2MGS{$valSI} = $MGS;
	}
	
	close OM;
	#die "$outD\n";
	if ($MGStax eq ""){
		print "No MGStax given, exiting without Abundance matrix calc..\n";
		exit (0);
	} 
	return \%Si2MGS;
}


sub transferSI2MGS{
	my ($hr) = @_;
	my %Si2MGS = %{$hr};
	foreach my $valSI (keys %Si2MGS){
		my $MGS = $Si2MGS{$valSI};
		#print "$valSI  $MGS\n";
		#get tax transferred..
		$specIfullTax{$MGS} = $specIfullTax{$valSI} unless (exists( $specIfullTax{$MGS} ));
		#actually replace SI with this MGS
		foreach my $COG (keys %{$SpecIgenes{$valSI}}){
			#transfer Q2S, delete old entry for it...
			foreach my $gen (@{$SpecIgenes{$valSI}{$COG}}){
				$Q2S{$gen}{$MGS} = $Q2S{$gen}{$valSI};
				delete $Q2S{$gen}{$valSI};
			}
			my @adds= @{$SpecIgenes{$valSI}{$COG}};
			if (exists($SpecIgenes{$MGS}{$COG})){ 
				push(@adds,@{$SpecIgenes{$MGS}{$COG}});
				@{$SpecIgenes{$MGS}{$COG}} = do { my %seen; grep { !$seen{$_}++ } @adds };
			} else { @{$SpecIgenes{$MGS}{$COG}} = @adds;}
			delete $SpecIgenes{$valSI}{$COG};
		}
		delete $SpecIgenes{$valSI};
	}
}


sub rebase($){ #calculates the profile for each SI based on the 40 marker genes
	my ($hr1)=@_; 
	my %specIs = %{$hr1};
	my $medianC=0;
	foreach my $sp (keys %specIs){
	#last;
		
		my @tar;my $MGn=0;
		foreach my $gid (@{$specIs{$sp}}){
			$MGn++;
			#now get genes from matrix
			#print "$gid\n";
			die "FMG entry missing: $gid\n" unless (exists($FMGmatrix{ $gid }));
			if ($medianC){
				for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){push(@{$tar[$j]}, ${$FMGmatrix{ $gid }}[$j]);}
			} else {
				for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$tar[$j] +=  ${$FMGmatrix{ $gid }}[$j];}
			}
		}
		if ($medianC){
			my @tar2;for (my $j=0;$j<scalar(@tar);$j++){$tar2[$j] = median(@{$tar[$j]});} $specIprofiles{$sp} = \@tar2;
		} else {
			for (my $j=0;$j<scalar(@tar);$j++){$tar[$j] /= $MGn;}$specIprofiles{$sp} = \@tar;
		}
		# if ($MGn>2);
	}
}

sub sanityCheckCorr(){
	#final sanity check, that marker genes are correlating
	my $wrongGene=0; my $corrGene=0;
	my %specIGcorrs;my %specIcnts; my %specIs;
	my @genes = keys %gene2specI;
	foreach my $k (@genes){
		my $g2si0  = (keys %{$gene2specI{$k}})[0];
		$specIcnts{$g2si0}++;push(@{$specIs{$g2si0}},$k);
	}
	
	#rebase..
	rebase(\%specIs);
	
	
	#actual correlation check
	foreach my $k (keys %gene2specI){
		#$specIGset{$gene2specI{$k}}
		my @sIS = keys %{$gene2specI{$k}};
		foreach my $sI (@sIS){ #in case several SIs have been assigned to gene..
			next unless ($specIcnts{$sI} > 2);
			next unless (exists($specIprofiles{$sI}));
			#if (nonZero($specIprofiles{$sI}) < 3){ next;}
			#print "@{$FMGmatrix{ $k }}\n";
			my $corr = correlation($specIprofiles{$sI},$FMGmatrix{ $k }) ;
			#print $corr." ";
			if ($corr < $globalCorrThreshold){
				$wrongGene++;rm4geneList($sI,$k);$specIcnts{$sI}--;
				next;
			}
			$corrGene++;
			push (@{$specIGcorrs{$sI}},$corr);
		}
	}
	print "correct(corr): $corrGene, rem(corr): $wrongGene. \n";
	
	undef %specIs;
	foreach my $k (keys %gene2specI){push(@{$specIs{   join(",",keys %{$gene2specI{$k}})   }},$k);}
	rebase(\%specIs);
}

sub passBlast($ $ $){
	my $spl = shift;my $reqID = shift; my $inM8 = shift;
	my $lengthGood = 0;
	if (@{$spl} < 12 || !defined($spl->[13])|| !defined($spl->[12])){ #prob serious error.. delete file
		print "not enough entries:\n@{$spl}\n$inM8\n";
		return -2;
	}
	$lengthGood=1 if ($spl->[3] >= ($spl->[13]*0.1) && $spl->[3] >= ($spl->[12]*0.7)  ); #subject query
	if ($spl->[2] < ($reqID*0.985) || !$lengthGood)  {return 0;}
	return 1;
}
sub lambdaBl($ $ $){
	my ($tar,$DB, $taxblastf) = @_;
	my $cmd="";
	if ($BlastCores > 20){$BlastCores = 20;}
	if (!-d $DB.".lambda/"  ) {
		print "Building LAMBDA index anew (may take up to an hour)..\n";
		my $cmdIdx = "$lambdaIdxBin -p blastn -t ".int($BlastCores)." -d $DB";
		if (system($cmdIdx)){die ("Lamdba ref DB build failed\n$cmdIdx\n");}
	}
	$cmd .= "$lambdaBin -t $BlastCores -id 93  -nm 100 -p blastn -e 1e-40 -q $tar -oc \"std qlen slen\" -i $DB.lambda -o $taxblastf\n";
	#die $cmd."\n";
	systemW $cmd if (!-e $taxblastf || $reblast);
}


sub readMotuTax($){
	my ($inF) = @_;
	my %gene2motu;
	my %motu2tax;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2motu{$spl[0]} = $spl[8];
		$motu2tax{$spl[8]} = join(";",@spl[1,2,3,4,5,6,7]) if (!exists $motu2tax{$spl[8]} );
	}
	close I;
	return (\%gene2motu,\%motu2tax);
}

sub readGene2mlinkage($){
	my ($inF) = @_;
	my %gene2LG;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2LG{$spl[0]} = $spl[2];
	}
	close I;
	return (\%gene2LG);
}
sub read_speci_tax($){
	my ($inF) = @_;
	my %gene2LG;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2LG{$spl[0]} = $spl[2];
	}
	close I;
	return (\%gene2LG);
}

sub readNCBI($){
	my ($inF) = @_;
	my %gene2LG;
	open I,"<$inF";
	while (my $l = <I>){
		chomp $l; my @spl=split/\t/,$l;
		$gene2LG{$spl[0]} = $spl[2];
	}
	close I;
	return (\%gene2LG);
}



sub specImatrix($ $){
	my ($oF,$hr) = @_;
	my %sTax = %{$hr};
	#print "@{$sTax{specI_v2_Cluster34}}\n";
	
	#create background count of SpecI genes not assigned
	my @bkgrnd; my @dblCh;
	foreach my $gid (keys %FMGmatrix){
		next if ($gid eq "header");
		for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$dblCh[$j] +=  ${$FMGmatrix{ $gid }}[$j]/40;}
		next if (exists ($gene2specI{$gid}));
		for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$bkgrnd[$j] +=  ${$FMGmatrix{ $gid }}[$j] / 40;}
	}
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	
	open Ox,">$oF" or die "Can't open out mat $oF\n";
	print Ox "SpecI\t".join ("\t",@{$FMGmatrix{ header }})."\n";
	#print O "SUM\t\t".join ("\t",@dblCh)."\n";
	print Ox "?\t".join ("\t",@bkgrnd)."\n";
	foreach my $si (sort keys %specIprofiles){
		print Ox "$si\t". join ("\t",@{$specIprofiles{ $si }})."\n";
	}
	close Ox;
	my $oFx = $oF; $oFx =~ s/\.[^\.]*$//;
	
	#and print SI tax for later reference..
	open Ot,">$oFx.tax" or die "Can;t open tax file $oFx.tax\n";
	foreach my $si (sort keys %specIprofiles){
		if (!exists($sTax{$si})){
			if (exists($MGSlist{$si})){
			} else {die "doesnt exist: $si\n";	}
		}
		print Ot "$si\t".join ("\t",@{$sTax{$si}})."\n";
	}
	close Ot;
	
	#
	#
	#calculating higher level abundance matrix
	my @taxLs = ("superkingdom","phylum","class","order","family","genus","species");
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	for (my $t=0;$t<@taxLs;$t++){
		#sum up to hi lvl
		my %thisMap;
		foreach my $si (sort keys %specIprofiles){
			
			if (!exists($sTax{$si})){
				if (exists($MGSlist{$si})){
				} else {
					die "doesnt exist: $si\n";
				}
			}
			print "ERR: $si    @{$sTax{$si}}\n" if (@{$sTax{$si}} <= $t);
			my $clvl = join (";",@{$sTax{$si}}[0 .. $t]);
			
			#if ($clvl eq "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia albertii"){
				#print "$si\n@{$specIprofiles{ $si }}\n";
			#}
			
			if (exists($thisMap{$clvl})){
				for (my $kl=0;$kl<scalar(@{$specIprofiles{ $si }});$kl++){
					${$thisMap{$clvl}}[$kl] += ${$specIprofiles{ $si }}[$kl];
				}
			} else {
				$thisMap{$clvl} = [@{$specIprofiles{ $si }}];
			}
		}
		open Ot,">$oFx.$taxLs[$t]" or die "Can't openn $oFx.$taxLs[$t]\n" ;
		print Ot "$taxLs[$t]\t".join ("\t",@{$FMGmatrix{ header }})."\n";
		print Ot "?\t".join ("\t",@bkgrnd)."\n";
		foreach my $kk (sort keys %thisMap){
			chomp $kk;
			print Ot "$kk\t".join("\t",@{$thisMap{$kk}}) . "\n";
		}
		close Ot;
	}
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";

}
sub nonZero($){
	my ($ar) = @_;
	my @a = @{$ar};
	my $nc =0 ;
	foreach my $x (@a){ $nc ++ if ($x > 0);}
	return $nc;
}

sub writeSpecItax{
	my ($hr1,$file) = @_;#\$specIid,"$SpecID/specI.tax2");
	my %sNTID = %{$hr1};
	open O,">$file" or die "can't open $file\n";
	foreach my $k (keys %sNTID){
		print O $k."\t".join("\t",@{$sNTID{$k}})."\n";
	}
	close O;
	die;
}

sub fixGTDBtax{
	my @spl=@_;
	if ($spl[1] eq "Firmicutes"){
		my $replTerm="Firmicutes_A";
		if ($spl[2] eq "Tissierellia" || $spl[2] eq "Erysipelotrichia"){$spl[2] = "Bacilli_A";$replTerm="Firmicutes_I";}
		if ($spl[2] eq "Bacilli"){$replTerm="Firmicutes_I";$spl[2] = "Bacilli_A";}
		if ($spl[2] eq "Peptococcia"){$replTerm="Firmicutes_B";}
		if ($spl[2] eq "Negativicutes"){$replTerm="Firmicutes_C";}
		if ($spl[2] eq "Bacilli"){$replTerm="Firmicutes_I";}
		$spl[1]=$replTerm;
	}
	if ($spl[1] eq "Proteobacteria" && $spl[2] eq "Deltaproteobacteria"){
		$spl[1] = "Desulfobacterota_A"; $spl[2] = "Desulfovibrionia";
	}
	if ($spl[1] eq "Bacteroidetes"){$spl[1] = "Bacteroidota";}
	if ($spl[5] eq "Methanobrevibacter"){$spl[5]="Methanobrevibacter_A";}
	if ($spl[3] eq "Campylobacterales"){$spl[1] = "Campylobacterota"; $spl[2] = "Campylobacteria";}
	if ($spl[4] eq "Pasteurellaceae"){$spl[3]="Enterobacterales";}
	if ($spl[4] eq "Succinivibrionaceae"){$spl[3]="Enterobacterales";}
	if ($spl[4] eq "Eggerthellaceae"){$spl[3] = "Coriobacteriales";}
	if ($spl[4] eq "Bifidobacteriaceae"){$spl[3]="Actinomycetales";$spl[1] = "Actinobacteriota";}
	if ($spl[3] eq "Actinomycetales"){$spl[1] = "Actinobacteriota";}
	if ($spl[2] eq "Actinobacteria"){$spl[2] = "Actinomycetia";}
	
	if ($spl[2] eq "Actinobacteria"){$spl[1] = "Actinobacteriota";}
	if ($spl[2] eq "Coriobacteriia"){$spl[1] = "Actinobacteriota";}
	if ($spl[1] eq "Actinobacteria"){$spl[1] = "Actinobacteriota";}

	if ($spl[5] eq "Aeromonas"){$spl[3]="Enterobacterales";}
	if ($spl[5] eq "Anaerococcus"){$spl[4]="Helcococcaceae";$spl[1]="Firmicutes_A";$spl[2]="Clostridia";}
	if ($spl[5] eq "Barnesiella"){$spl[4]="Barnesiellaceae";}
	if ($spl[5] eq "Blautia"){$spl[3]="Lachnospirales";}
	if ($spl[5] eq "Bradyrhizobium"){$spl[4]="Xanthobacteraceae";}
	if ($spl[5] eq "Caecibacter"){$spl[3]="Veillonellales";$spl[4]="Megasphaeraceae";}
	if ($spl[5] eq "Chryseobacterium"){$spl[2]="Bacteroidia";$spl[4]="Weeksellaceae";}
	if ($spl[5] eq "Corynebacterium"){$spl[3]="Mycobacteriales";$spl[4]="Mycobacteriaceae";}
	if ($spl[5] eq "Dialister"){$spl[4]="Dialisteraceae";}
	if ($spl[5] eq "Dorea"){$spl[4]="Lachnospirales";}
	if ($spl[5] eq "Dysgonomonas"){$spl[4]="Dysgonomonadaceae";}
	if ($spl[5] eq "Eubacterium"){$spl[3]="Eubacteriales";}
	if ($spl[5] eq "Eubacterium_E"){$spl[4]="Lachnospiraceae";}
	if ($spl[5] eq "Leuconostoc"){$spl[4]="Lactobacillaceae";}
	if ($spl[5] eq "Megasphaera"){$spl[4]="Megasphaeraceae";}
	if ($spl[5] eq "Mogibacterium"){$spl[3]="Peptostreptococcales";$spl[4]="Anaerovoracaceae";}
	if ($spl[5] eq "Pantoea"){$spl[4]="Enterobacteriaceae";}
	if ($spl[5] eq "Parabacteroides"){$spl[4]="Tannerellaceae";}
	if ($spl[5] eq "Paraprevotella"){$spl[4]="Bacteroidaceae";}
	if ($spl[5] eq "Peptostreptococcus"){$spl[3]="Peptostreptococcales";}
	if ($spl[5] eq "Prevotella"){$spl[4]="Bacteroidaceae";}
	if ($spl[5] eq "Pseudoflavonifractor"){$spl[3]="Oscillospirales";$spl[4]="Oscillospiraceae";}
	if ($spl[5] eq "Ruminococcus"){$spl[3]="Oscillospirales";}
	if ($spl[5] eq "Serratia"){$spl[4]="Enterobacteriaceae";}
	if ($spl[5] eq "Staphylococcus"){$spl[3]="Staphylococcales";}
	if ($spl[5] eq "Butyricicoccus"){$spl[3]="Oscillospirales";$spl[4]="Butyricicoccaceae";}
	if ($spl[5] eq "Campylobacter"){$spl[1]="Campylobacterota";$spl[2]="Campylobacteria";}
	if ($spl[5] eq "Lachnoclostridium"){$spl[3]="Lachnospirales";}
	if ($spl[5] eq "Paenibacillus"){$spl[3]="Paenibacillales";}
	if ($spl[5] eq "Slackia"){$spl[3]="Coriobacteriales";}	return @spl;
}

#parsing of specI taxonomy
sub createAreadSpecItax{
	my ($hr1,$file,$GTDBtax) = @_;#\$specIid,"$SpecID/specI.tax2");
	my %sNTID = %{$hr1}; my $tFileHits=0;
	my @taxIds = keys %sNTID;
	my %ret; my $cnt = -1;
	if (-e $file){
		open I,"<$file";
		my $cnt = 0;
		while (my $l = <I>){
			$cnt++;
			chomp $l; my @spl = split /\t/,$l;
			#if ($cnt == 0){}
			my $id = shift @spl;
			while (@spl < 7){
				push(@spl,"?");
			}
			if ($doGTDBtax){ #GTDB has Fimicutes_A etc
				@spl = fixGTDBtax(@spl);
			}
			for (my $i=0;$i<@spl;$i++){if ($spl[$i] eq "" || $spl[$i] eq " "){$spl[$i] ="?";}}
			$ret{$id} = \@spl;
			$tFileHits++;
			$cnt++;
		}
		close I;
		print "Read $cnt precompiled entries\n";
	}
	
	if ($tFileHits == scalar @taxIds){return(\%ret);}
	
	
	#print "$ret{specI_v2_Cluster1309}\n";
	my %GTDB;
	if (-e $GTDBtax){
		print "Reading GTDB taxids\n";
		open I,"<$GTDBtax" or die $!;
		while (<I>){
			chomp; my @spl = split /\t/;
			my $gt1 = $spl[0]; $gt1 =~ s/[dpcofgs]__//g;
			my @gt = split /;/,$gt1;
			#die "@gt\n";
			$GTDB{$spl[1]}=\@gt;
		}
		close I;
	}
	
	my @taxids; my @spids;
	open O,">>$file"; my $GTDBhits=0;
	foreach my $k (keys %sNTID){
		next if (exists($ret{$k}));
		#die "$sNTID{$k} $k\n";
		if (exists($GTDB{$sNTID{$k}})){
	#take tax directly from GTDB
			print O "$k\t".join("\t",@{$GTDB{$sNTID{$k}}})."\n";
			$GTDBhits++;
		} else {
			push (@taxids,$sNTID{$k});
			push (@spids,$k);
		}
	}
	close O;
	print "\nFound $GTDBhits GTDB hits\n";
	
	
	#die @taxids."\n";
	#die "python /g/bork3/home/hildebra/dev/python/get_ranks.py ".join(" ",@taxids);
	if (@taxids>0){
		print "Detecting ".@taxids." new taxids\n";
		my $pyTaxId = getProgPaths("taxid2tax_scr");
		my $cmd = "$pyTaxId ".join(" ",@taxids);
		my $tret= `$cmd`;
		my @newT = split /\n/,$tret;
		if (@newT > 0){
			open O,">>$file";
			for (my $i=0;$i<@newT;$i++){
				my @spl = split /\t/,$newT[$i]; shift @spl;
				print O "$spids[$i]\t".join("\t",@spl)."\n";
			}
			close O;
		}
	}
	#die "@{$ret{specI_v2_Cluster34}}\n";
	return (\%ret);
}

sub rm4geneList($ $){
	my ($k,$sg) = @_;
	#die "ASDA";
	delete $gene2specI{$sg}{$k};
	if (keys(%{$gene2specI{$sg}}) == 0 ){delete $gene2specI{$sg};
	} else {#more complicated, remove from comma list
		#my @spl = keys %{$gene2specI{$sg}};
		#print "@spl\n";
		#@spl = grep (!/$k/i, @spl);
		#die "TODO: @spl  XX $k";
		
	}
	
	return unless (exists($SpecIgenes2{$k}));
	#print "A";
	#$c is unknown..
	my $entrFound=0;
	foreach my $c (keys %{$SpecIgenes2{$k}}){
		#print "L";
		if ($SpecIgenes2{$k}{$c} eq $sg){
			delete $SpecIgenes2{$k}{$c}; $entrFound=1;last;
		}
	}
	#die "not deleted from SpecIgenes2 $k $sg\n" unless ($entrFound==1);
	#return $ret;
}


#check which mean abundance profile single MG have, and how multi MGs correlate to this 
#then selects best correlating MG to be "the one" that just fits
sub getCorrs{
	my $dblAssi=0;
	my $dblA=0;my $singlA=0;my $singlMultA=0;my $skippedSIs=0; my $newAssigns=0;
	foreach my $k(keys %SpecIgenes){ #this is specI
		my @tarGenes; #matrix vector summed over single MGs
		my $gcnt =0; my %selV;
		next if (exists($MGSlist{$k})); #this is an MGS and doesn't have Q2S etc and should be static in any case..
		foreach my $c(keys %{$SpecIgenes{$k}}){ #this is COG
			$selV{$c} =0;
			if (@{$SpecIgenes{$k}{$c}}>1){next;}#multi copy, dont use this gene
			my $gid = ${$SpecIgenes{$k}{$c}}[0];
			if (scalar (keys (%{$Q2S{$gid}})) > 1 ){$dblA++;next;}
			$singlA++; 
			#if ($gcnt == 0){@tarGenes = @{$FMGmatrix{ $gid }};
			#} else {
			for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$tarGenes[$j] += ${$FMGmatrix{ $gid }}[$j];}
			#} 
			$gcnt++;				
			#$gene2specI{$gid} = $k;   #can still be wrong assignment.. 
			$dblAssi++ if (add2geneList($k,$c,$gid));
			$selV{$c} =1;
			#if (exists($gene2specI{$gid})){
			#	$gene2specI{$gid} .= ",".$k; $dblAssi++;#print "should not happen\n";#$gene2specI{$spl[$i]}   $k\n";
			#} else {$gene2specI{$gid} = $k; } 
		}
		#print "T2: $gcnt   " . scalar(keys(%selV))." D" if ($k eq "TEC2");
		
		
		if ($gcnt < 3 && scalar(keys %{$SpecIgenes{$k}}) >= 3){ #no use, calc everything together up
			#These genes will later be checked again for the correlation to the grand mean
			foreach my $c(keys %{$SpecIgenes{$k}}){
			#last;
				my @spl = @{$SpecIgenes{$k}{$c}};
				if (@spl > 1){next;}
				my @subTarG; my $sgcnt=0;
				foreach my $gid (@spl){#single copy, use this gene
					if (scalar (keys (%{$Q2S{$gid}})) > 1){next;}#don't want multi species assignments
					#if ($sgcnt == 0){@subTarG = @{$FMGmatrix{ $gid }}; 
					#} else {
						for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$subTarG[$j] += ${$FMGmatrix{ $gid }}[$j];}
					#}
					$sgcnt++;
				} 
				if ($sgcnt>0){#norm by number added genes
					$singlMultA++;
					$gcnt++ ;
					for (my $j=0;$j<@subTarG;$j++){$tarGenes[$j] += $subTarG[$j] / $sgcnt;}
				}
			}

		}
		if ($gcnt ==0){$skippedSIs++;next;}	
		if (nonZero(\@tarGenes) < 3){ next;}#print "XX" ;
		#doesn't need norm, since we do spearman correlation
		#but now corr and see which genes just fit best of the multi choices..
		#print "@tarGenes\n";
		foreach my $c(keys %{$SpecIgenes{$k}}){
			if (!$selV{$c}){ #$SpecIgenes{$k}{$c} =~ m/,/){#only look at multi assigned genes
				#die "$SpecIgenes{$k}{$c}\n";
				my @spl = @{$SpecIgenes{$k}{$c}};
				my @subCors; my $max=0;
				foreach my $sg (@spl){
					die "$sg doesn't exist in gene list \n" if (!exists($FMGmatrix{ $sg }));
					
					#my $corr = calculate_spearman_correlation(\@tarGenes,$FMGmatrix{ $sg }) ;
					#print "@{$FMGmatrix{ $sg }}\n";
					my $corr = correlation(\@tarGenes,$FMGmatrix{ $sg }) ;
					if ($corr > $max){$max = $corr;}
					push (@subCors,  $corr  );
				}
				if ($max < $globalCorrThreshold){next;}
				$newAssigns++;
				for (my $i=0;$i<@spl;$i++){
					unless ($subCors[$i]>$max-0.03){next;}
					my $sg = $spl[$i];
					#this assignment is what I need, now I know this gene is blocked for assignment to other SpecI's
					$dblAssi++ if (add2geneList($k,$c,$sg));
					for (my $j=0;$j<@tarGenes;$j++){$tarGenes[$j] += ${$FMGmatrix{ $sg }}[$j] ;}
				}
			}
		}

		#norm vector
		for (my $j=0;$j<@tarGenes;$j++){$tarGenes[$j] /= $gcnt;}
		#and save the final specI profile..
		$specIprofiles{$k} = \@tarGenes;
		#print @tarGenes ."XX\n";
	}
	print "double assignment $dblAssi; assigned: $newAssigns   Stats in Run: mult. spec. $dblA $singlA $singlMultA $skippedSIs\n";
}

#just create the abundances for MGS related FMG's
sub createProfileMGS{
	foreach my $k(keys %SpecIgenes){ #this is specI
		my @tarGenes; #matrix vector summed over single MGs
		my $gcnt = scalar keys %{$SpecIgenes{$k}}; 
		next unless (exists($MGSlist{$k})); #this is an MGS and doesn't have Q2S etc and should be static in any case..
		#print "\n$k ";
		foreach my $c(keys %{$SpecIgenes{$k}}){
			#die "$SpecIgenes{$k}{$c}\n";
			my @spl = @{$SpecIgenes{$k}{$c}};
			$gcnt++;
			#print "@spl ";
			for (my $i=0;$i<@spl;$i++){
				my $sg = $spl[$i];
				#this assignment is what I need, now I know this gene is blocked for assignment to other SpecI's
				#$dblAssi++ if (add2geneList($k,$c,$sg));
				#for (my $j=0;$j<@tarGenes;$j++){$tarGenes[$j] += ${$FMGmatrix{ $sg }}[$j] ;}
				for (my $j=0;$j<scalar(@{$FMGmatrix{ $sg }});$j++){$tarGenes[$j] += ${$FMGmatrix{ $sg }}[$j];}
			}
		
		}

		for (my $j=0;$j<@tarGenes;$j++){$tarGenes[$j] /= $gcnt;}
		#add MGS profile
		$specIprofiles{$k} = \@tarGenes;
	}
}

sub disentangleMultiAssigns(){ 
#some genes are assigned to more than one specI; these need to be checked for better corr
	foreach my $k (keys %gene2specI){
		my @spl = keys %{$gene2specI{$k}};
		if (@spl > 1 ){
			my @refAB = @{$FMGmatrix{ $k }};
			
			my $max = 0;  my $maxID=0;
			my @subCorr; my @subIDs;
			foreach my $sI (@spl){
				my $corr = correlation(\@refAB,$specIprofiles{ $sI }) ;
				$max = $corr if ($corr>$max);
				my $gID = $gen2SIscore{$k}{$sI};
				$maxID = $gID if ($gID > $maxID);
				push(@subCorr,$corr);
				push(@subIDs,$gID);
			}
			my $nzc = nonZero(\@refAB);
			my $subhits=0;my $best= "";
			for (my $i=0;$i< @spl ; $i++){
				if ($subCorr[$i] >= ($max-0.01)){#hit
					$subhits++;
					$best = $spl[$i];
				}
			}
			if ($subhits > 1 || $nzc < 4){
				$subhits=0;
				for (my $i=0;$i< @spl ; $i++){
					if ($subIDs[$i] >= ($maxID-0.01)){#hit
						$subhits++;
						$best = $spl[$i];
					}
				}
			}
			if ($subhits>1){
				print"too many subhits $k : $subhits\n@spl\n@subCorr\n@subIDs\n";
			}
			if ($best ne "" ){
				rmFromList($k,$best);
				$gene2specI{$k} = {$best=>undef};
			}

		}
	}
}


sub tree4FMGs{
	my $btout = "$GCd/FMG/specIphylo/";
	system "mkdir -p $btout" unless (-d $btout);
	my $hr = readFasta("$GCd/FMG/COG*.faa"); my %FAA = %{$hr};
	$hr = readFasta("$GCd/FMG/COG*.fna"); my %FNA = %{$hr};
	my $SaSe = "|";
	open ON,">$btout/all.fna"; open OA,">$btout/all.faa"; 
	my %catT;
	foreach my $SI (keys %SpecIgenes){
		my @CGs = keys %{$SpecIgenes{$SI}};
		foreach my $cg (@CGs){
			my $gID = ${$SpecIgenes{$SI}{$cg}}[0];
			next unless (exists($gene2specI{$gID})); #just make sure this specI is also in the latest FMG counting...
			$gID =~ m/^([^,]+)/; $gID = $1;
			die "@{$SpecIgenes{$SI}{$cg}}\n" unless (exists($FNA{$gID}));
			print ON ">$SI$SaSe$cg\n$FNA{$gID}\n";
			print OA ">$SI$SaSe$cg\n$FAA{$gID}\n";
			push(@{$catT{$cg}},"$SI$SaSe$cg");
		}
	}
	close ON; close OA;
	open OC,">$btout/all.cats";
	foreach my $cg (keys %catT){
		print OC join("\t",@{$catT{$cg}})."\n";
	}
	close OC;
	print "Creating phylogeny for found specI's//\n";
	my $cmd= "$bts  -aa  $btout/all.faa -smplSep '\\$SaSe' -cats $btout/all.cats -outD $btout -runIQtree 1 -runFastTree 0 -runRaxMLng 1 -cores $BlastCores  -AAtree 1 -bootstrap 000 -NTfiltCount 300 -NTfilt 0.1 -NTfiltPerGene 0.5 -minOverlapMSA 2 -MSAprogram 2 -AutoModel 0 -iqFast 1 \n";
	systemW "$cmd\n";
	my $QSBoptHR = emptyQsubOpt(1,"");
	#my ($dep,$qcmd) = qsubSystem($btout."treeCmd.sh",$cmd,$BlastCores,"1G","FMGtree","","",1,[],$QSBoptHR);
}

sub rmFromList{ #just removes from SpecIgenes2 list..
	my ($g,$s) = @_;
	foreach my $c (keys %{$SpecIgenes2{$s}}){
		if (defined $SpecIgenes2{$s}{$c} && $SpecIgenes2{$s}{$c} eq $g){
			delete $SpecIgenes2{$s}{$c} ;
			last;
		}
	}
}

sub add2geneList($ $ $){ #assign a gene ($sg) to a speci($k), and its COG ($c)
	my ($k,$c,$sg) = @_;
	my $ret=0;
	if (exists($gene2specI{$sg})){
		 $ret=1;#print "should not happen\n";#$gene2specI{$spl[$i]}   $k\n";
	} 
	$gene2specI{$sg}{$k} = undef;
	$SpecIgenes2{$k}{$c}=$sg;#set mark to block this MG in this specI...
	return $ret;
}





