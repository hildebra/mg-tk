#!/usr/bin/env perl
#annotates all marker genes in gene catalog


use warnings;
use strict;

use Mods::GenoMetaAss qw( readClstrRev systemW median readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt qsubSystemJobAlive);
use Mods::geneCat qw(calculate_spearman_correlation read_matrix correlation);
use Mods::TamocFunc qw( readTabbed3);
use Mods::FuncTools qw(passBlast lambdaBl);


use Mods::IO_Tamoc_progs qw(getProgPaths);
use List::Util;
use Getopt::Long qw( GetOptions );

sub reformatGTDBtax;
sub submitJobs; sub missingLCAs;

my $LCAbin = getProgPaths("LCA");

#0.1: 7.12.23: added lambda3 support
my $version = 0.1;

my $GCd ="";#$ARGV[0]."/";
my $cores = 1;#$ARGV[1];
my $useGTDBmg = "FMG"; #GTDB or FMG
my $tmpD = "";
my $qsubDir = "";
my $totMem = "50"; #30G should be enough..

#options to pipeline..
GetOptions(
	"GCd=s"      => \$GCd,
	"tmp=s"      => \$tmpD,
	"cores=i"    => \$cores,
	"MGset=s"    => \$useGTDBmg,#GTDB or FMG
);


#submission system setup
$qsubDir = $GCd."/LOGandSUB/MarkerGenes/";
system "mkdir -p $qsubDir" unless (-d $qsubDir);
my $doSubmit = 1;my $submSys = "";

my $QSBoptHR = emptyQsubOpt($doSubmit,"",$submSys);#,"bash"
$QSBoptHR->{qsubDir} = $qsubDir;



#prep MGS related paths..
my $speciesLink = "specI_lnks"; my $speciesCutoff = "specI_cutoff";
my $speciesGTDB = "specI_GTDB"; my $speciesDir = "specIPath";my $subsetFile = "FMG.subset.cats";
my $MGtag = "FMG";
#GTDB marker genes instead of FMG?? 
die "-MGset option has to be \"GTDB\" or \"FMG\"\n" unless ($useGTDBmg eq "GTDB" || $useGTDBmg eq "FMG");
if ($useGTDBmg eq "GTDB"){ 
	$speciesLink = "GTDB_lnks"; $speciesCutoff = "GTDB_cutoff"; $subsetFile = "GTDBmg.subset.cats";
	$speciesGTDB = "GTDB_GTDB"; $speciesDir = "GTDBPath"; $MGtag = "GTDBmg";
}
my $inSImap = getProgPaths($speciesLink);
my $GTDBspecI = getProgPaths($speciesGTDB);
my $SpecID=getProgPaths($speciesDir);#directoy with all 40 SpecI marker genes
my %FMGcutoffs = %{readTabbed3(getProgPaths($speciesCutoff,0),1)};

#make the distinction here, as not all Marker genes might be present in a certain experiment (e.g. no Archaea)
my %FMGkeys = %{readTabbed3("$GCd/$subsetFile",1)};





my $taxPerGene = "$SpecID/$MGtag.tax";
reformatGTDBtax($taxPerGene, $inSImap,$GTDBspecI);
#die;

#set up file structure
system "mkdir -p $tmpD";# or die "MG_LCA.pl:: can't create tmpdir $tmpD\n";


my $MGdir = "$GCd/$MGtag/";
my $redo = 0; my $subm=1;
#	lambdaBl($tar,$DB,$taxblastf);
my $blPerJob = 40;my $xtrLab = "";

my $jobrounds=0;
while (missingLCAs() ){
	submitJobs();
	$jobrounds++;
	if ($jobrounds>2){
		die "Something really wrong, more than 2 jobrounds..";
	}
}

system "rm -r $tmpD";
print "Done with $MGtag LCAs\n\n"; 

system " cat $MGdir/*.LCA | cut -f2,3,4,5,6,7,8 | sort | uniq -c | sort -r -n -k1 > $MGdir/$MGtag.freqs.txt";

exit(0);




sub missingLCAs{
	
	my $miss =0;
	foreach my $COG (keys %FMGkeys){#(@catsPre){
		my $ifna = "$MGdir/$COG.fna";
		my $CogTaxF = $ifna; $CogTaxF =~ s/\.fna$/\.LCA/;
		$miss++ unless (-e $CogTaxF);
	}
	
	print "Missing $miss of " . scalar(keys %FMGkeys) . " MG LCAs\n";
	return $miss;
}



sub submitJobs{
	my $COGcnt = 0; my @Deps;
	my $cmd = ""; my $collectJobs=0;
	#print "Start subm\n";
	my $avx2Constr = getProgPaths("avx2_constraint",0);
	my @preCons = @{$QSBoptHR->{constraint}};push(@{$QSBoptHR->{constraint}}, $avx2Constr);
	foreach my $COG (keys %FMGkeys){#(@catsPre){
		my $ifna = "$MGdir/$COG.fna";
		die "Could not find input fna: $ifna\n" unless (-e $ifna);
		system ("ln -s $MGdir/$COG.fna $MGdir/$COG.fa") unless (-e "$MGdir/$COG.fa");
		$ifna = "$MGdir/$COG.fa";
		#die $ifna." $SpecID/$COG${xtrLab}.fna $tmpD/$COG${xtrLab}.tmp.m8\n";
		my $m8file = "$tmpD/$COG${xtrLab}.tmp.m8";
		system "rm $m8file" if ($redo);
		my $CogTaxF = $ifna; $CogTaxF =~ s/\.fn?a$/\.LCA/;
		unless (-e $CogTaxF || -e $m8file){
			#print "XX\n";
			$cmd .= "\n\n#At $COG\n";
			$cmd .= lambdaBl($ifna,"$SpecID/$COG${xtrLab}.fna",$m8file,$cores,!$subm)."\n" unless (-e $CogTaxF || -e $m8file); #.rep
			$cmd .= "$LCAbin  -i $m8file -r $taxPerGene -o $CogTaxF  -LCAfrac 0.8  -cover 0.9 -minAlignLen 70 -id $FMGcutoffs{$COG},90,80,60,50,30,0;\n" unless (-e $CogTaxF);
			$collectJobs++;
		}
		$COGcnt++;
		#die $cmd."\n$m8file\n";
		
		if ( ($COGcnt % $blPerJob) == 0 && $collectJobs>= $blPerJob){
			my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
			my ($dep1,$qcmd1) = qsubSystem($qsubDir."MG_LCA.$COGcnt.sh",$cmd,$cores,int($totMem/$cores)."G","MLCA$COGcnt","","",1,[],$QSBoptHR);
			$QSBoptHR->{tmpSpace} =$tmpSHDD;
			push (@Deps,$dep1);
			#die;
			$cmd = "";
			$collectJobs=0;
			print "Submitted $COG blast ($COGcnt / ". scalar(keys(%FMGcutoffs)) . ")\n";
		}
	}
	#submit remainder..
	if ($cmd ne ""){
		print "SUB\n\n";
		my $tmpSHDD = $QSBoptHR->{tmpSpace};	$QSBoptHR->{tmpSpace} = "0"; 
		my ($dep1,$qcmd1) = qsubSystem($qsubDir."MG_LCA.$COGcnt.sh",$cmd,$cores,int($totMem/$cores)."G","MLCA$COGcnt","","",1,[],$QSBoptHR) ;
		$QSBoptHR->{tmpSpace} =$tmpSHDD;
		push (@Deps,$dep1);
		#print $cmd."\n";
		$cmd = "";
	}
	@{$QSBoptHR->{constraint}} = @preCons;
	#print "Submitted $COGcnt / ". scalar(keys(%FMGcutoffs)) . " \n";

	print "Waiting for MG LCA's..\n";
	qsubSystemJobAlive( \@Deps,$QSBoptHR ); 
}




sub reformatGTDBtax{
	my ($taxPerGene, $inSImap,$GTDBspecI ) = @_;
	return if (-e $taxPerGene);
	my %sp2id = %{readTabbed3($inSImap,1)};
	open O,">$taxPerGene";
	open I, "<$GTDBspecI";
	while (my $line = <I>){
		chomp $line; my @spl = split /\t/,$line;
		my $spec = shift @spl; shift @spl;
		if (!exists $sp2id{$spec}){
			die "MG_LCA::Can't find $spec\n";
		}
		print O $sp2id{$spec} . "\t" . join (";",@spl) . "\n";
	}
	close I;
	close O;
	print "Wrote $taxPerGene\n";
	#die;
}
