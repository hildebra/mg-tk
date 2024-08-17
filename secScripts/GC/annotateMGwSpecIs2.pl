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

use Mods::GenoMetaAss qw( readClstrRev systemW median mean readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt );
use Mods::geneCat qw(calculate_spearman_correlation read_matrix correlation);
use Mods::FuncTools qw(passBlast lambdaBl);
use Mods::TamocFunc qw( readTabbed3 readTable);
use Mods::math qw(nonZero);


use Mods::IO_Tamoc_progs qw(getProgPaths);
use List::Util;
use Getopt::Long qw( GetOptions );

#.11: reverted some of the drastic gene filtering
#.12 adapted for Bin_SB folder potentially being different..
my $version = 0.12;

sub readMotuTax; sub fixGTDBtax;
sub readGene2mlinkage; sub readSpecIids;
sub readNCBI;sub read_speci_tax;
sub createProfileMGS;
sub MGSassign; sub transferSI2MGS;
sub readMGS;
sub rebase0;

#sub calculate_spearman_correlation;
#sub read_matrix; 

sub getCorrs;
sub add2geneList; sub rm4geneList;
sub sanityCheckCorr;
sub specImatrix;
sub rmFromList;
#sub createAreadSpecItax; 
sub writeSpecItax;
sub disentangleMultiAssigns;

#my $SpecID="/g/bork3/home/hildebra/DB/MarkerG/specI/"; my $freeze11=1;
my $freeze11=0;
my $doGTDBtax = 1; #Firmicutes_A etc
#progenomes.specIv2_2
my $globalCorrThreshold = 0.6; # determines cutoff, when still to accept correlating genes into specI
my $reblast=0;#do blast again?

my $rarBin = getProgPaths("rare");#"/g/bork5/hildebra/dev/C++/rare/rare";
my $samBin = getProgPaths("samtools");#"/g/bork5/hildebra/bin/samtools-1.2/samtools";
my $bts = getProgPaths("buildTree_scr");



if (@ARGV == 0){
	die "Not enough input args: use ./annotateMGwMotus.pl [path to GC] [# Cores]\n";
}

my $GCd ="";#$ARGV[0]."/";
my $mode = "specI";
my $BlastCores = 1;#$ARGV[1];
my $MGSfile = ""; #MGS annotations; will take precedence of blast specI annotations..
#$MGSfile = $ARGV[2] if (@ARGV > 2);
my $MGStax = ""; #MGS taxonomy (fitting to MGSfile).. only if given will try to calc abundance table
#$MGStax = $ARGV[3] if (@ARGV > 3);
my $useGTDBmg = "FMG"; #GTDB or FMG
my $tmpD = "";
my $minGenes = 10; #that many MGs are required, to include a species..
my $outD = "";
my $hr1; #general purpose hash ref pointer..




#options to pipeline..
GetOptions(
	"GCd=s"      => \$GCd,
	"outD=s"     => \$outD,
	"tmp=s"      => \$tmpD,
	"MGS=s"      => \$MGSfile,
	"cores=i"    => \$BlastCores,
	"MGStax=s"   => \$MGStax,
	"MGset=s"    => \$useGTDBmg,#GTDB or FMG
	"minGenes=i" => \$minGenes,
);
die "Needs option -GCd $GCd\n" if ($GCd eq "");


my $speciesLink = "specI_lnks"; my $speciesCutoff = "specI_cutoff";
my $speciesGTDB = "specI_GTDB"; my $speciesDir = "specIPath";
my $MGterm = "FMG";
#GTDB marker genes instead of FMG?? 
die "-MGset option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");
if ($useGTDBmg eq "GTDB"){ 
	$MGterm = "GTDBmg";
	$speciesLink = "GTDB_lnks"; $speciesCutoff = "GTDB_cutoff";
	$speciesGTDB = "GTDB_GTDB"; $speciesDir = "GTDBPath";
}
my %FMGcutoffs = %{readTabbed3(getProgPaths($speciesCutoff,0),1)};
my $MGdir = "$GCd/$MGterm/";
$outD = "$GCd/Anno/Tax/$MGterm/" if ($outD eq "");
if ($MGSfile ne ""){$outD = "$GCd/Anno/Tax/${MGterm}_MGS/";}
my $inSImap = getProgPaths($speciesLink);
my $GTDBspecI = getProgPaths($speciesGTDB);
my $SpecID=getProgPaths($speciesDir);#directoy with all 40 SpecI marker genes


print "-------------------------------------------------------------------------\nannotateMG script v $version\n-------------------------------------------------------------------------\n";
print "Using gene cat in $GCd\nMGS: $MGSfile\nMGStax: $MGStax\ncores: $BlastCores MGset: $useGTDBmg minGenes: $minGenes\n";


#ystem "mkdir -p $MGdir/tax" unless (-d "$MGdir/tax");
system "mkdir -p $outD" unless (-d "$outD");

#die "Writing to $outD, using $MGdir\n";

#die $outD;
my $motuDir = "";#"/g/bork3/home/hildebra/DB/MarkerG/mOTU";

#annotate against DB using lambda
#my $hr1 = readTable($GTDBspecI,"\t",";" );
#my %specItax = %{$hr1};


#first check if taxpergene exists (can be used in LCA algo directly)
if (0 ){ #not needed: there should be *.LCA file now!
	my $taxPerGene = "$SpecID/specI.pergene.tax";
	my $prepPG = getProgPaths("progenomes_prep_scr");
	my $cmd = "$prepPG $SpecID $inSImap $GTDBspecI $taxPerGene\n";
	system $cmd if (!-e $taxPerGene);
}

my ($hrID, $hrTAX) = readSpecIids($inSImap);
my %specIid = %{$hrID};#my %specItax = %{$hrTAX};

#specItax is essentially only used to register what key is present..
#my $specIfullTaxHR = createAreadSpecItax(\%specItax,"$SpecID/specI.tax3",$GTDBspecI); #specI.tax3 doesn't need to exist any longer..
#
#tax per specI - new way relying on precomputed LCAs
my %specIfullTax = %{readTable($GTDBspecI,"\t" )};# %{$specIfullTaxHR};
my $xtrLab= "";$xtrLab= ".rep" if ($freeze11);

#assign each COG separately
system "mkdir -p $MGdir" unless (-d $MGdir);

my %specItaxname; my %SpecIgenes; 
my %COGDBLass; my %COGass; my %gene2COG;
my %gen2SIscore;

my %Q2S; #assignment of GC genes to specI id (poss. several ids)

my %COG2gen; 
#my %FMGlist;
my $COGgenes=0;
open IC,"<$GCd/${MGterm}.subset.cats" or die "Can't open $GCd/${MGterm}.subset.cats\n";
while (<IC>){
	chomp;	my @spl  = split /\t/;
	#$cats{$spl[0]} = $spl[2];
	my @genes = split(/,/,$spl[2]);
	my $curCOG = $spl[0];
	$COG2gen{$curCOG} = \@genes;
	foreach (@genes){$gene2COG{$_} = $curCOG; $COGgenes++;}
}
close IC;
print "Read $COGgenes COG 2 MGs genes\n";




#MGS related containers
my %Gene2MGS; my %speci2MGS; 
my %MGSlist; #$MGSlist{$curMGS} = 1;
my %speci2MGScnt; #counts how many different genes were used in total

readMGS($MGSfile); #only used if in MGS mode..
#undef %FMGlist; #no longer needed from here on
#die;
my $allOK=1;

#new routine: will read .LCA to get links..
print "Reading LCA's\n";
foreach my $COG (keys %COG2gen){#(@catsPre){
	my @genes = @{$COG2gen{$COG}};
	#my $finiN=0;my %finiM=();
	#print "Reading $MGdir/$COG.LCA\n";
	$hr1 = readTable("$MGdir/$COG.LCA","\t",";",1);
	my %tax = %{$hr1};
	foreach my $gen (keys %tax){
		my $speci = $tax{$gen}; #complete tax string.. well should be ok
		#$gene2COG{$gen} = $COG;
		if (!exists($specIfullTax{$speci})){
			my @tmp = split /;/,$speci;
			$specIfullTax{$speci} = \@tmp;
		}
		
		#from here old code..
		#just link the specI to MGS..
		if (exists($Gene2MGS{$gen})){
			my $curMGS = $Gene2MGS{$gen};
			$speci2MGS{$curMGS}{$speci}++;
			#if ($finiN == 0){$speci2MGScnt{$curMGS} ++ ; $finiN=1;} $finiM{$speci}=1;
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
}

#die;

if ($allOK==0){die"One or more Blast files were not ok.. restart procedure\n";}

#write out SpecI -> MGS assignments (and tax)
my $SI2MGShr = MGSassign();
#transfer gene assignments to MGS, transfer SI tax -> MGS tax
transferSI2MGS($SI2MGShr);
my %SI2MGS = %{$SI2MGShr};


#die;
#print "Done initial blast\n$MGdir\n";
#foreach (sort {$COGass{$a} cmp $COGass{$b}} keys %COGass){print "$_ $FMGcutoffs{$_}:  $COGass{$_}($COGDBLass{$_})\n";}
print "Reading MG marker gene matrix..\n";
$hr1 = read_matrix("$GCd/Matrix.$MGterm.mat"); #here I need the coverage matrix..
#my $hr1 = read_matrix("$GCd/Mat.cov.FMG.mat"); #here I need the coverage matrix..
#my $hr1 = read_matrix("$GCd/Mat.med.FMG.mat"); #here I need the coverage matrix..
my %FMGmatrix= %{$hr1};

my %gene2specI; my %specIprofiles; 
my %SpecIgenes2; #collects list of genes that could be associated to specI
#sort out best multi hit by correlation analysis
getCorrs();#\%FMGmatrix,\%SpecIgenes);


rebase0();#create first specIprofile..
#print "SIZE:: " . scalar(keys %gene2specI) . "\n";
disentangleMultiAssigns();

#my $T2cnt=0;foreach my $k (keys %gene2specI){$T2cnt++ if ($gene2specI{$k} eq "TEC2");}print "T2cnt = $T2cnt\n";

#check that all corrs here check out well..
#sanityCheckCorr();

#print "DEBUG "; foreach my $k (keys %specIprofiles){ print "$k ; \n" if ($k =~ m/,/);}

print "resolving remainder genes..\n";

undef %Q2S; my $xtraEntry=0;
#add genes though correlations.. ??? actually uesless..
foreach my $COG (keys %COG2gen){#(@catsPre){
	#last;
	my %specIcnt;
	my @genes = @{$COG2gen{$COG}};
	my $reqID = $FMGcutoffs{$COG};
	#secondary scanning of assignments...
#LCA based code..
	#print "Reading $MGdir/$COG.LCA\n";
	my %tax = %{readTable("$MGdir/$COG.LCA","\t",";",1)};
	
	foreach my $gid (keys %tax){
		next if (exists($Gene2MGS{$gid}));
		my $speci = $tax{$gid}; #complete tax string.. well should be ok
		$speci = $SI2MGS{$speci} if (exists($SI2MGS{$speci}));
		#3 this MG has already been assigned in the high confidence initial assignments
		if (exists($SpecIgenes2{$speci}{$COG} ) ){next;} #|| exists($Gene2MGS{$gid})
		
		#4 correlate to species core, to make sure the gene kind of makes sense..
		if (exists($specIprofiles{$speci})){
			my $corr = correlation($specIprofiles{$speci},$FMGmatrix{ $gid }) ;
			next if ($corr < $globalCorrThreshold);
			#print $corr."\t";
		} else { next;}

		if (!exists($Q2S{$gid}{$speci} )){
			$Q2S{$gid}{$speci} = scalar(keys(%{$Q2S{$gid}})); #print "m";
		}

		$specIcnt{$speci}++;
		$xtraEntry++;
		#and block gene slot
		add2geneList($speci,$COG,$gid); # == $gene2specI{$gid}{$speci} = undef;

	}
	#print "$COG: Found ".keys(%Q2S)." assignments (".@genes.")\n";
	
}

#disentangleMultiAssigns();

#my $k = "4070197";my @spl2 = keys %{$gene2specI{$k}};print "BEF = @spl2\n";print "\n\n@spl2\n" ;

print "Entries not in MGS: $xtraEntry\n";


#write specI assignments for markerG
open O,">$MGdir/gene2specI.txt";
foreach my $k (keys %gene2specI){
	if (!exists($gene2COG{$k})){
		print "No COG assignment $k!  ";
	} else {
		print O "$k\t". join(",",keys %{$gene2specI{$k}}) . "\t$gene2COG{$k}\n";
	}
}
close O;

#add MGS profiles extra to the -> no longer needed..
#createProfileMGS();
#create abundance profile
rebase0();#create final specIprofile..
specImatrix("$outD/specI.mat",\%specIfullTax);


print "Finished SpecI annotations & matrix: $outD\n";
exit(0);




#####################################################
#####################################################
#####################################################
#####################################################



sub readMGS{
	my ($MGSfile) = @_;
	return if ($MGSfile eq "");
	print "Reading reference MGS: $MGSfile\nTax: $MGStax\n";
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
			#@spl = fixGTDBtax(@spl) if ($doGTDBtax);
			$tmpTax{$id} = \@spl;
		}
		close IT;
	}
	#print"T::@{$tmpTax{MGS0287}}\n";
	my $taxF=0; my $taxN=0;
	open IM,"<$MGSfile" or die "Can't open MGS $MGSfile\n";
	while (<IM>){ 
		chomp; my @spl = split /\t/; my @spl2= split /,/,$spl[1]; 
		my $curMGS = $spl[0]; my $genCnt=0;
		foreach my $gen (@spl2){
			#chomp $gen;
			next unless (exists($gene2COG{$gen}));
			#print "$gen " if ($curMGS eq "MB2bin186");
			$Gene2MGS{$gen} = $curMGS;
			$MGSlist{$curMGS} = 1;
			#just insert in specI related objects.. just pretend it's great
			#die "$gen $curMGS\n" if ($gen <5);
			push(@{$SpecIgenes{$curMGS}{$gene2COG{$gen}}},$gen);
			$gen2SIscore{$gen}{$curMGS}=200;
			#$genCnt++;
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
		#print "0 genes found for MGS $curMGS\n";
	}
	close IM;
	my $meanS=0; my @sizes; my $morethan1=0; my $only1=0;
	foreach my $MGS (keys %SpecIgenes){
		#print $MGS." ";
		my $curS =0;
		foreach my $cat (keys %{$SpecIgenes{$MGS}}){
			 my $lcurS = @{$SpecIgenes{$MGS}{$cat}};
			 $curS += $lcurS;
			 if ($lcurS > 1){$morethan1 ++ ;} else {$only1++;}
		}
		$meanS += $curS;
		push(@sizes, $curS);
	}
	@sizes = sort @sizes;
	print "Found ". scalar(keys %MGSlist)." MGS with $taxF/".($taxF+$taxN) ." taxonomies. Mean MGs/MGS: " . $meanS /($taxF+$taxN) . "; median: " . median(@sizes) .". $morethan1/".($only1+$morethan1)." with >1 copy.\n";
#die "@sizes\n";
	
	#print"T::@{$specIfullTax{MGS0287}}\n";

}


sub MGSassign{
	#$speci2MGS{}
	my %Si2MGS ;
	return \%Si2MGS if ($MGSfile eq "");
	my $logfile = "$outD/MGS2speci.txt";
	open OM,">$logfile";
	my $notAssigned=0; my $assignedMGS =0; my @maxAssi;
	foreach my $MGS (sort {$speci2MGS{$b} <=> $speci2MGS{$a}} keys %speci2MGS){
		print OM "$MGS";
		my @SIsAdd;
		my $valSI="";my $valStr=0;
		my $totalAssi =0;
		my $maxassFrac = 0;
		foreach my $si (keys %{$speci2MGS{$MGS}}){$totalAssi += $speci2MGS{$MGS}{$si}; }
		
		foreach my $si (sort {$speci2MGS{$MGS}{$b} <=> $speci2MGS{$MGS}{$a}} keys %{$speci2MGS{$MGS}}){
			 push(@SIsAdd,"$si:$speci2MGS{$MGS}{$si}");
			 if ($valSI eq "" && $speci2MGS{$MGS}{$si}>5){
				$valSI = $si;$valStr =$speci2MGS{$MGS}{$si}; 
			 } elsif ($valStr < 0.5*$speci2MGS{$MGS}{$si}){ #just too ambigous..
				$valSI = "";
			 }
			my $assFrac= $speci2MGS{$MGS}{$si}/$totalAssi;
			$maxassFrac = $assFrac if ($assFrac >= $maxassFrac);
		}
		push(@maxAssi,$maxassFrac);
		if ($valSI ne ""){
			#print for file the link
			print OM "\t".join(";",@{$specIfullTax{$valSI}});
			$assignedMGS++;
		} else {
			$notAssigned++;
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
	
	#some more assignment stats..
	my $chimera=0;
	foreach (@maxAssi){
		$chimera ++ if ($_ <0.6);
	}
	my $totMGS = $notAssigned+$assignedMGS;
	print "Matched ${assignedMGS} of ". ($totMGS)  ." MGS to taxa/specIs: $logfile\n";
	#some stats on how uniform LCA was..
	print "Mean max assignments MGS: " . mean(@maxAssi) . " , median: " .median(@maxAssi) ."; pot. chimeric (species level and above): $chimera \n";
	
	
	
	return \%Si2MGS;
}


sub transferSI2MGS{
	my ($hr) = @_;
	my %Si2MGS = %{$hr};
	print "Comparing MGS to SpecIs..\n" if (scalar(keys %Si2MGS)>0);
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
	my $medianC=0; my $skippedSpecies=0;
	foreach my $sp (keys %specIs){
	#last;
		#print "$sp- ";
		my @tar;my $MGn=0;
		foreach my $gid (@{$specIs{$sp}}){
			$MGn++;
			#now get genes from matrix
			#print "$gid\n";
			die "${MGterm} entry missing: $gid\n" unless (exists($FMGmatrix{ $gid }));
			if ($medianC){
				for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){push(@{$tar[$j]}, ${$FMGmatrix{ $gid }}[$j]);}
			} else {
				for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){$tar[$j] +=  ${$FMGmatrix{ $gid }}[$j];}
			}
		}
		#print "$sp: " . scalar(@tar) . " $MGn \n";
		if (0 && $MGn < $minGenes){
			$skippedSpecies++;
			next;
		}
		
		if ($medianC){
			my @tar2;for (my $j=0;$j<scalar(@tar);$j++){$tar2[$j] = median(@{$tar[$j]});} 
			$specIprofiles{$sp} = \@tar2;
		} else {
			for (my $j=0;$j<scalar(@tar);$j++){$tar[$j] /= $MGn;}
			$specIprofiles{$sp} = \@tar;
		}
		# if ($MGn>2);
	}
	print "Rebbase: skipped $skippedSpecies species/MGS\n";
}

sub rebase0{
	my %specIs;
	foreach my $k (keys %gene2specI){
		my @spl = sort keys %{$gene2specI{$k}};
		die "gene not in gene2specI $k\n" if (!exists($gene2specI{$k}) || @spl == 0);
		my $g2si0  = $spl[0];
		#die "$k   $g2si0\n"; = 4124843   Bacteria;Firmicutes;Bacilli;E..
		push(@{$specIs{$g2si0}},$k);
	}
	rebase(\%specIs);
}

sub sanityCheckCorr(){
	#sanity check, that marker genes are correlating
	return;
	my $wrongGene=0; my $corrGene=0;
	my %specIGcorrs;my %specIcnts; my %specIs;
	foreach my $k (keys %gene2specI){
		my $g2si0  = (keys %{$gene2specI{$k}})[0];
		$specIcnts{$g2si0}++;
		push(@{$specIs{$g2si0}},$k);
	}
	
	#rebase..
	rebase(\%specIs);
	
	
	#actual correlation check
	foreach my $k (keys %gene2specI){
		#$specIGset{$gene2specI{$k}}
		next if (exists($Gene2MGS{$k}));
		foreach my $sI (keys %{$gene2specI{$k}}){ #in case several SIs have been assigned to gene..
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
	print "Creating specI matrix..\n";
	my %sTax = %{$hr};
	#print "@{$sTax{specI_v2_Cluster34}}\n";
	
	my %specIcnts; #just use for histgram
	foreach my $k (keys %gene2specI){$specIcnts{ (keys %{$gene2specI{$k}})[0] }++;}
	
	my $rmSpecs=0;my %delSIs;
	foreach my $si (sort keys %specIprofiles){
		if (exists($specIcnts{$si}) && $specIcnts{$si} < $minGenes){
			delete $specIprofiles{$si};
			$delSIs{$si} = 1;
			$rmSpecs++;
		}
	}
	print "Removed $rmSpecs specI's from final matrix due to having <$minGenes marker genes\n";
	
	#create background count of SpecI genes not assigned
	my %bkgrnd; my @dblCh;
	foreach my $gid (keys %FMGmatrix){
		next if ($gid eq "header");
		for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){       $dblCh[$j] +=  ${$FMGmatrix{ $gid }}[$j] ;} #TODO:40
		if (exists ($gene2specI{$gid})  ){
			if ( exists( $delSIs{  (keys %{$gene2specI{$gid}})[0] }  )    ){
				delete $gene2specI{$gid};# if (exists($gene2specI{$gid}));
				#delete $gene2specI{$k} if (exists($gene2specI{$k}));
			} else {
				next;
			}
		}
		if (exists ($Gene2MGS{$gid}) ){next;}
		for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){      $bkgrnd{$gene2COG{$gid}}{$j} +=  ${$FMGmatrix{ $gid }}[$j] ;}
	}
	#split up backgrnd by gene id/sample
	my @bkgrnd1;
	for (my $j=0;$j<scalar(@dblCh);$j++){ 
		my @meanVal ; my $ccnt=0;
		foreach my $cog (keys %bkgrnd){
			if (exists($bkgrnd{$cog}{$j}) && $bkgrnd{$cog}{$j} > 0){
				$ccnt++; push(@meanVal, $bkgrnd{$cog}{$j});
			}
		}
		
		#$meanVal /= $ccnt if ($ccnt > 0);
		$bkgrnd1[$j] = mean(@meanVal);
	}
	
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	
	open Ox,">$oF" or die "Can't open out mat $oF\n";
	print Ox "SpecI\t".join ("\t",@{$FMGmatrix{ header }})."\n";
	#print O "SUM\t\t".join ("\t",@dblCh)."\n";
	print Ox "?\t".join ("\t",@bkgrnd1)."\n";
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
		print Ot "?\t".join ("\t",@bkgrnd1)."\n";
		foreach my $kk (sort keys %thisMap){
			chomp $kk;
			print Ot "$kk\t".join("\t",@{$thisMap{$kk}}) . "\n";
		}
		close Ot;
	}
	#print "@{$specIprofiles{specI_v2_Cluster34}}\n";
	
	
	
	#some stats on created abundances.. 
	%specIcnts = ();
	foreach my $k (keys %gene2specI){$specIcnts{ (keys %{$gene2specI{$k}})[0] }++;}
	my %histo;
	for my $k (sort {$specIcnts{$a} <=> $specIcnts{$b}} keys %specIcnts) {
		$histo{int $specIcnts{$k}/10}++;
	   # print "$k $specIcnts{$k} $NTax{$specItax{$k}}\n" ;#if ($specIcnts{$k}>=40);   # bbb c aaaa
	}
	foreach (sort{$a <=> $b} (keys %histo)){
		print "$_\t$histo{$_}\n";
	}


}

sub writeSpecItax{
	my ($hr1,$file) = @_;#\$specIid,"$SpecID/specI.tax2");
	print "Writing newly created  SpeciI tax to $file\n";
	my %sNTID = %{$hr1};
	open O,">$file" or die "can't open $file\n";
	foreach my $k (keys %sNTID){
		print O $k."\t".join("\t",@{$sNTID{$k}})."\n";
	}
	close O;
	die;
}

sub readSpecIids($){
	my ($inSImap ) = @_;
	#specItax is essentially only used to register what key is present..
	my %specIid;my %specItax;
	open I,"<$inSImap" or die "Can't open in specI map $inSImap\n";
	#specI v3 (progenomes2)
	#while (<I>){next if (m/^#/);chomp; my @xx = split /\t/;$xx[1] =~s/,//g; $specIid{$xx[1]} = $xx[0];$xx[1]=~m/^(\d+)\./; $specItax{$xx[0]}=$1;}
	#specI v4 (progenomes3)
	while (<I>){next if (m/^#/);chomp; 
		my @xx = split /\t/;
		foreach my $yy (split /;/,$xx[1]){
		$yy =~s/,//g; $specIid{$yy} = $xx[0];
		$xx[1]=~m/^(\d+)\./; $specItax{$xx[0]}=$1;
		}
	}
	close I;
	return (\%specIid, \%specItax);
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
	my $belowGeneIncl = 0;
	print "Correlations of MGs..\n";
	foreach my $k (keys %SpecIgenes){ #this is specI
		print "$k ; " if ($k =~ m/,/);
		my @tarGenes; #matrix vector summed over single MGs
		my $gcnt =0; my %selV;
		
		#commented, as gene needs to be stored in $SpecIgenes2
		#next if (exists($MGSlist{$k})); #this is an MGS and doesn't have Q2S etc and should be static in any case..
		
		foreach my $cog(keys %{$SpecIgenes{$k}}){ #this is COG
			$selV{$cog} =0;
			if (@{$SpecIgenes{$k}{$cog}}>1){next;}#multi copy, dont use this gene
			my $gid = ${$SpecIgenes{$k}{$cog}}[0];
			if (scalar (keys (%{$Q2S{$gid}})) > 1 ){$dblA++;next;}
			$singlA++; 
			#print " $gid ";
			#if ($gcnt == 0){@tarGenes = @{$FMGmatrix{ $gid }};
			#} else {
			for (my $j=0;$j<scalar(@{$FMGmatrix{ $gid }});$j++){
				$tarGenes[$j] += ${$FMGmatrix{ $gid }}[$j];
			}
			#} 
			$gcnt++;				
			$dblAssi++ if (add2geneList($k,$cog,$gid));
			$selV{$cog} =1;
		}
		#print "T2: $gcnt   " . scalar(keys(%selV))." D" if ($k eq "TEC2");
		
		
		if (0 && $gcnt < 3 && scalar(keys %{$SpecIgenes{$k}}) >= 3){ #no use, calc everything together up
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
		#if ($gcnt < $minGenes){$belowGeneIncl++;next;}
		
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
	print "double assignment $dblAssi; assigned: $newAssigns   Stats in Run: mult. spec. $dblA $singlA $singlMultA $skippedSIs $belowGeneIncl\n";
	
	
}

#just create the abundances for MGS related FMG's
sub createProfileMGS{
	die "no longer used!!";
	foreach my $k(keys %SpecIgenes){ #this is specI
		next if (exists($specIprofiles{$k}));
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
	print "MGS profiles created\n";
}

sub disentangleMultiAssigns(){ 
	#some genes are assigned to more than one specI; these need to be checked for better corr
	my $resolvedDbl=0; my $problems =0 ;
	print "Disentangle genes double assignments..\n";
	foreach my $k (keys %gene2specI){
		 
		my @spl = sort keys %{$gene2specI{$k}};
		if (@spl <= 1 ){next;}
		$problems++;
		my $assign2MGS=0;
		
		
		die "Can't find gene $k in FMG matrix!\n" unless (exists($FMGmatrix{ $k }));
		my @refAB = @{$FMGmatrix{ $k }};
		
		my $max = 0;  my $maxID=0;
		my @subCorr; #my @subIDs;
		my @spl2;
		foreach my $sI (@spl){
			#die "Ref profile not found $sI ($k)\n" ;
			next unless (exists($specIprofiles{ $sI }));
			my $corr = correlation(\@refAB,$specIprofiles{ $sI }) ;
			my $gID = 1;
			if (exists($MGSlist{$sI})){
				$corr += 1;# = $gen2SIscore{$k}{$sI};
			}
			$max = $corr if ($corr>$max);
			#$maxID = $gID if ($gID > $maxID);
			push(@subCorr,$corr);
			push(@spl2,$sI);
			#push(@subIDs,$gID);
		}
		my $nzc = nonZero(\@refAB);
		my $subhits=0;my $best= "";
		for (my $i=0;$i< @spl2 ; $i++){
			if ($subCorr[$i] >= ($max-0.01) ){#hit
				$subhits++;
				$best = $spl2[$i];
				#$assign2MGS=1 if (exists($Gene2MGS{$k}));
			}
		}
		if ($subhits > 1 || $nzc < 4){
			$subhits=0;
			for (my $i=0;$i< @spl2 ; $i++){
				if ($subCorr[$i] >= ($max)){#hit
					$subhits++;
					$best = $spl2[$i];
				}
			}
		}
		if ($subhits>1){
			print"too many subhits $k : $subhits\n@spl2\n@subCorr\n";#@subIDs\n";
		}
		
		#ok has a prefered candidate.. just assign to this, unassign from others..
		if ($best ne "" ){
			#rmFromList($k,$best);
			#$gene2specI{$k} = {$best=>1};
			
			
			foreach (@spl2){
				if ($_ eq $best){
					add2geneList($_,$gene2COG{$k},$k);
				} else {
					rmFromList($_,$gene2COG{$k},$k);
				}
				
			}
			add2geneList($best,$gene2COG{$k},$k);
			$resolvedDbl++;
		}
	}
	print "Resovled $resolvedDbl genes asigned to MGS + specI (N=$problems problematic)\n";
}


sub tree4FMGs{
	my $btout = "$GCd/${MGterm}/specIphylo/";
	system "mkdir -p $btout" unless (-d $btout);
	my $hr = readFasta("$GCd/${MGterm}/*.faa"); my %FAA = %{$hr};
	$hr = readFasta("$GCd/${MGterm}/*.fna"); my %FNA = %{$hr};
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

sub rmFromList($ $ $){ #just removes from SpecIgenes2 list.. reverses "add2geneList"
	my ($k,$c, $sg) = @_;
	delete ( $SpecIgenes2{$k}{$c} );
	#print "RM $k, $c, $sg\n$SpecIgenes2{$k}{$c}\n";
	delete $gene2specI{$sg}{$k};# = 0;
	#foreach my $c (keys %{$SpecIgenes2{$sg}}){
	#	if (defined $SpecIgenes2{$sg}{$c} && $SpecIgenes2{$sg}{$c} eq $g){
	#		delete $SpecIgenes2{$sg}{$c} ;
	#		last;
	#	}
	#}
}

sub add2geneList($ $ $){ #assign a gene ($sg) to a speci($k), and its COG ($c)
	my ($k,$c,$sg) = @_;
	my $ret=0;
	if (exists($gene2specI{$sg})){
		 $ret=1;#print "should not happen\n";#$gene2specI{$spl[$i]}   $k\n";
	} 
	$gene2specI{$sg}{$k} = 1;
	$SpecIgenes2{$k}{$c}=$sg;#set mark to block this MG in this specI...
	return $ret;
}






