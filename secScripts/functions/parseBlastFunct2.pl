#!/usr/bin/env perl
#parseBlastFunct.pl /g/bork5/hildebra/data/Ach_proteins/res/DiaAssignment.txt [DB] [mode]
#mode == 0: just parse & report sum stat; 1=report per gene primary DB; 2=report per gene all levels ==2 check for all files present
use warnings;
use strict;
use FileHandle;
use Data::Dumper;

use Mods::FuncTools qw(readGene2COG);
use Mods::TamocFunc qw(sortgzblast readTabbed readTabbed2 readTabbed3 uniq);
use Mods::GenoMetaAss qw(gzipwrite gzipopen convertNT2AA);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Getopt::Long qw( GetOptions );


sub main;
sub readGene2COG; sub readNogKingdom; sub createNOGKgd;
sub readCOGdef; sub readTCDBdef;
sub readGene2KO;sub readKeggTax;
sub check_files; sub remove_file;
sub writeAllTable;
sub readMohTax;
sub combineBlasts;
sub bestBlHit;
sub help;
sub fixCardMeta;
sub readCzySubs;
sub eggMap_interpret;
my $emBin = getProgPaths("emapper");




my $blInf = "";#$ARGV[0];
my $mode = 0;#ARGV[3]
my $DBmode = "NOG";#$ARGV[1];
my $quCovFrac = 0; #how much of the subject (DB) needs to be covered?
my $noHardCatCheck = 0; #select for the hit with KO assignment rather than the real best hit (w/o KO assignment)

#$DBmode = uc $DBmode;
my $minBLE= 1e-7;#ARGV[2]
my $minScore=0;#min required bit score
my $minAlLen = 30; #my $minPidGlobal = 40;
my $tmpD="";#$ARGV[6]
my $DButil = "";#$ARGV[5];
my $lengthF = "";#$ARGV[4];
my $reportEggMapp=0; #do eggNOG mapper?
my $ncore = 4;
my $writeSumTbls = 1; #report all the higher level stats?
my $percID = 20;
my $Card_leniet = 1; #take more CARD hits: eval^(1/$Card_leniet)
my $calcGL = 1;
my $noTax=0;
my $KOfromNOG=0;
my $checkTaxNog = ""; #reads in NOG assignments of genes and determines tax from these (used for AB production, to check it's not bacterial gene)
my $writeFastaOut = 0; #extracts the domains that are matching to a hit
my $useBacNOG = 0; #use bacNOGs instead of LUCA level COGs

die "no input args!\n" if (@ARGV == 0 );
GetOptions(
	"help|?" => \&help,
	"i=s"      => \$blInf,
	"DB=s"      => \$DBmode,
	"mode=i"      => \$mode,#0=normal, 1=normal and print per gene anno, 2=print extended per gene annotation, no summary,3=file check 4=remove outputfile
	"eval=s"      => \$minBLE,
	"minBitScore=f" => \$minScore,
	"minAlignLen=i"      => \$minAlLen,
	"minPercSbjCov=f"      => \$quCovFrac,
	"minFractQueryCov=f"      => \$quCovFrac,
	"tmp=s"      => \$tmpD,
	"DButil=s"	=> \$DButil,
	"LF=s"	=> \$lengthF,
	"KOfromNOG=i"	=> \$KOfromNOG,
	"lenientCardAssignments=i" => \$Card_leniet,
	"eggNOGmap=i"	=> \$reportEggMapp,
	"summaryTbls=i" => \$writeSumTbls,
	"calcGeneLengthNorm=i" => \$calcGL,
	"singleSpecies=i" => \$noTax,
	"CPU=i" => \$ncore,
	"bacNOG=i" => \$useBacNOG,
	"reportDomains=i" => \$writeFastaOut,
	"NOGtaxChk=s" => \$checkTaxNog, 
	"percID=i" => \$percID, #percent id similiarity, from 0 - 100
) or die("Error in command line arguments\n");

die "Database mode requires -tmp -minAlignLen -LF\n" if ( ($mode != 4 && $mode != 3) && ( $tmpD eq "" || $lengthF eq "" ) );
#die "$writeSumTbls XX\n";
$reportEggMapp=0 if ($DBmode ne "NOG" );


my %DBlen;
if ($lengthF ne "" ){#read the length of DB proteins
	my @spl = split /,/,$lengthF;
	foreach my $lF2 (@spl){
		open I,"<$lF2" or die "length file $lF2 no present!\n";
		while (my $l = <I>){
			chomp $l;
			my @spl = split /\t/,$l;
			$DBlen{$spl[0]} = int $spl[1];
		}
		close I;
	}
	print "read length DB\n";
}

$blInf = sortgzblast($blInf,$tmpD) unless ($mode == 3 || $mode == 4);


#die "$blInf";

my $bl2dbF ;my $cogDefF ;my $NOGtaxf ;
my $KEGGlink ;my $KEGGtaxDb ; my $TCDBhir;
my $PATRVIRanno; my $VIRDBannp;
if ($DButil ne ""){
	$bl2dbF = "$DButil/NOG.members.tsv";
	$bl2dbF = "$DButil/bactNOG.members.tsv" if ($useBacNOG);
	$cogDefF = "$DButil/NOG.annotations.tsv";
	$NOGtaxf = "$DButil/all_species_data.txt";
	$KEGGlink = "$DButil/genes_ko.list";
	$KEGGtaxDb = "$DButil/kegg.tax.list";
	$TCDBhir = "$DButil/TCDBhir.txt";
	$PATRVIRanno = "$DButil/PATRIC_VF2.tab";
	$VIRDBannp = "$DButil/VF.tab";
}

$blInf =~ m/(.*)\/([^\/]+$)/;
my $inP = $1; my $inF = $2;
#die "$inP\t$inF\n";

my @blRes;
my @kgdOpts = qw (0 1 2 3);
my @kgdName = ("Bacteria","Archaea","Eukaryota","Fungi");
my @kgdNameShrt = ("BAC","ARC","EUK","FNG");
#which best blast decide algo to use..
my $scomode = 0; my $emode = 1;
my @normMethods = ("cnt");
if ($calcGL){push (@normMethods,"GLN");}
my $czyEuk =2; my $czyFungi=5;
my %czyTax;my %czySubs;
my @aminBLE = split /,/,$minBLE;
#my @aminPID = (40,45,50,55,60,65,70,75,80,85,90,95);
#@aminBLE = ("1e-9") x scalar(@aminPID);
my @aminPID = ($percID) x scalar(@aminBLE);

my %TCdef; my %PTVcat;my %PTVgene;my %PTVgene_def;
my %NOGkingd ; my %KEGGtax; my %read2NOGkgd;
my %COGdef ;my %g2COG ; my %c2CAT ; 
my %cardE; my %cardFunc; my %cardName; 
my %cardFull;my %cardAux;
my $tabCats = 0; my $cazyDB = 0; my $ACLdb = 0; my $KGMmode=0;

if ($checkTaxNog ne "" ){
	%read2NOGkgd = createNOGKgd($checkTaxNog);#0=BAC, 1= aRCH, 2=euks, 3=fung
	#die;
}
if ($mode == 3 || $mode == 4){ #scan for all reads finished
} elsif ( $DBmode eq "KGB" || $DBmode eq "KGE" || $DBmode eq "KGM"){
	print "Reading KEGG DBs..\n";
	my ($hr1) = readGene2KO($KEGGlink);
	%c2CAT = %{$hr1};
	#my @tmp = keys %c2CAT; die "$tmp[0] $tmp[1] $tmp[2]\n";
	@kgdOpts = qw (3); $tabCats = 2;
	if ($DBmode eq "KGB"){ @kgdName = ("","","","Bacteria");@kgdNameShrt = ("","","","BAC") ;}
	if ($DBmode eq "KGE"){ @kgdName = ("","","","Eukaryotes");@kgdNameShrt = ("","","","EUK") ;}
	if ($DBmode eq "KGM"){ @kgdName = ("Prokayotes","Eukaryotes","","Unclassified");@kgdNameShrt = ("BAC","EUK","","UNC") ; 
		@kgdOpts = qw (0 1 3); $KGMmode=1; 
	}
	$hr1 = readKeggTax($KEGGtaxDb);
	%KEGGtax = %{$hr1};
	system "ln -s $DButil/ko_descr/ko_KO.descr $inP/KEGG.descr" unless (-e "$inP/KEGG.descr");

} elsif ( $DBmode eq "NOG"){
	print "Reading NOG DBs..\n";
	my ($hr1,$hr2) = readGene2COG($bl2dbF);
	%g2COG = %{$hr1}; %c2CAT = %{$hr2};
	$hr1 = readCOGdef($cogDefF);
	%COGdef = %{$hr1}; $tabCats = 1;
	%NOGkingd = readNogKingdom($NOGtaxf);
	system "ln -s $DButil/NOG.anno $inP/NOG.anno" unless (-e "$inP/NOG.anno");
	system "ln -s $DButil/NOG.descr $inP/NOG.descr" unless (-e "$inP/NOG.descr");


} elsif ($DBmode eq "TCDB" || $DBmode eq "PTV"|| $DBmode eq "VDB" ){
	@kgdOpts = qw (3);
	@kgdName = ("","","","ALL");
	@kgdNameShrt = ("","","","ALL");
	$tabCats =6;
	if ($DBmode eq "TCDB"){
		$tabCats =6;
		my $hr = readTCDBdef($TCDBhir);
		system "ln -s $TCDBhir $inP/TCDB.anno" unless (-e "$inP/TCDB.anno");
		%TCdef = %{$hr};
		print "transporter DB\n";
	}
	if ($DBmode eq "VDB"){
		$tabCats = 8;
		print "VIRDB DB\n";
		my $hr = readTabbed3($VIRDBannp,4);
		%PTVcat = %{$hr};
		$hr = readTabbed3($VIRDBannp,1);
		%PTVgene = %{$hr};
		$hr = readTabbed3($VIRDBannp,2);
		%PTVgene_def = %{$hr};
		system "ln -s $VIRDBannp $inP/VirDB.anno" unless (-e "$inP/VirDB.anno");

	}
	if ($DBmode eq "PTV"){
		$tabCats=8;
		print "PATRIC virulence factors DB\n";
		my $hr = readTabbed3($PATRVIRanno,3);
		%PTVcat = %{$hr};
		$hr = readTabbed3($PATRVIRanno,1);
		%PTVgene = %{$hr};
		$hr = readTabbed3($PATRVIRanno,2);
		system "ln -s $PATRVIRanno $inP/PTV.anno" unless (-e "$inP/PTV.anno");

		%PTVgene_def = %{$hr};
		foreach my $k (keys %PTVgene){
			if ($PTVgene{$k} eq ""  ){
				if ($PTVgene_def{$k} ne ""){
					$PTVgene_def{$k} =~ m/^(\S+)/;
					$PTVgene{$k} = $1;
				} else {
					$PTVgene{$k} = "unknwn";
				}
			}
		}
		#die $PTVcat{'fig|984235.3.peg.1162|VBISalEnt284594_1162|'}."\n";
	}

} elsif ($DBmode eq "CZy"){
	@kgdOpts = qw (0 1 2 3 4 5);
	my $hr = readMohTax($DButil."/MohCzy.tax");
	%czyTax = %{$hr};
	$hr = readCzySubs("$DButil/cazy_substrate_info.txt");
	#die "ln -s $DButil/CAZyDB.07312019.fam-activities.txt $inP/Cazy.families";
	#system "ln -s $DButil/CAZyDB.07312019.fam-activities.txt $inP/Cazy.families" unless (-e "$inP/Cazy.families");
	system "ln -s $DButil/cazy_substrate_info.txt $inP/cazy_substrate_info.txt" unless (-e "$inP/cazy_substrate_info.txt");
	#system "ln -s $DButil/cazy_substrate_reformated.txt $inP/cazy_substrate_info2.txt" unless (-e "$inP/cazy_substrate_info2.txt");
	
	%czySubs = %{$hr};
	$czyEuk =2; $czyFungi=5;
	#2015
	#@kgdName = ("Bacteria","Archaea","Eukaryota","unclassified","CBM","Fungi");
	#2019
	@kgdName = ("Bacteria","Archaea","Eukaryota","NA","Viruses","Fungi");
	@kgdNameShrt = ("BAC","ARC","EUK","UNC","VIR","FNG");
	$cazyDB = 1;
	$tabCats=3;
	print "CAZy database\n";
} elsif ($DBmode eq "ACL"){
	@kgdOpts = qw (3);
	$czyEuk =2; $czyFungi=5;
	@kgdName = ("Plasmid","Prophage","Virus","unclassified");
	@kgdNameShrt = ("PLA","PRO","VIR","UNC");
	$ACLdb = 1;
	print "aclame database\n";
} elsif ($DBmode eq "PAB"){
	@kgdOpts = qw (3);
	@kgdName = ("","","","Fungi");
	@kgdNameShrt = ("","","","FNG");
	print "Fungi AntiBiotics production domains DB\n";
	$tabCats=7;
} elsif ($DBmode eq "ABRc"){
	@kgdOpts = qw (3);
	@kgdName = ("","","","All");
	@kgdNameShrt = ("","","","ALL");
	print "ABR Card database\n";
	if (0){#old style dependant on full length gene predictions
		my $hr = readTabbed($DButil."cardEval.txt");
		%cardE = %{$hr};
		$hr = readTabbed($DButil."cardClass.txt");
		%cardFunc = %{$hr};
		$hr = readTabbed($DButil."cardName.txt");
		%cardName = %{$hr};
	} else {
		my ($hr,$MaxDepth) = readTabbed2($DButil."card.parsed.f11.tab.map",0);
		system "ln -s $DButil/card.parsed.f11.tab.map $inP/Card.tax" unless (-e "$inP/Card.tax");
		$hr = fixCardMeta($hr);
#		($hr,$MaxDepth) = readTabbed2($DButil."card.parsed.f11.txt",0);
		%cardFull = %{$hr};
#		%cardAux = %{$hr};
		
		#print $cardFull{"GI:CAP12351.2"}."\n";
		#die $cardFull{"796891.PRJNA163057.ATYY01000007_gene147"}."\n";

	}
	$tabCats=5;
} else{
	@kgdOpts = qw (3);
	@kgdName = ("","","","All");
	@kgdNameShrt = ("","","","ALL");
	print "Non - NOG DB\n";
}

if ($noTax){
	@kgdOpts = qw (3);
	@kgdName = ("","","","All");
	@kgdNameShrt = ("","","","ALL");
	print "no Tax info being used\n";
}


#read DB file
my $readLim = 1000000; my $lcnt=0;my $lastQ=""; my $stopInMiddle=0; my $reportGeneCat=0;
my $mergeDiaHit =0;  my $noMerge=0;
my %COGhits; my %CAThits; my %GENEhits; my $O2; my %OEM;
my %rds2xtrct;
my %emFiles ; my %emFilesFinal;
my %FAO;

for (my $i=0; $i<@aminBLE ; $i++){$COGhits{$i} = {}; $CAThits{$i} = {}; $GENEhits{$i} = {};}

if ($mode == 0 || $mode==1 || $mode == 2){ #mode1 = write gene assignment, mode 2: write for each gene higher lvl cat
	#die "$blInf\n";
	$reportGeneCat = $mode if ($mode >= 1);
	#insert here KEGG BAC ? EUK fix
	if ($writeFastaOut){
		foreach my $x (@kgdOpts){ #open outstream for fasta output
			for (my $i=0; $i<@aminBLE ; $i++){
				my $pathXtra = "/CNT_".$aminBLE[$i]."_".$aminPID[$i]."/";
				my $outPath = $inP.$pathXtra;
				system "mkdir -p $outPath" unless (-d $outPath);
				open $FAO{$i}{$x},">$outPath/$DBmode.faa";
			}
		}
	}
	if ($KGMmode && !-e $blInf){
		#zcat dia.KGB.blast.gz dia.KGE.blast.gz | sort | gzip > dia.KGM.blast.gz
		my $KGBf = $blInf; my $KGEf = $blInf; 
		$KGBf =~ s/KGM/KGB/; $KGEf =~ s/KGM/KGE/;
		my $cmd = "zcat $KGBf $KGEf | sort | gzip > $blInf\n";
		system $cmd."\n";
	}
	
	my ($I,$OK) = gzipopen($blInf,"diamond output file [in]",1); 
	my $OK2;
	($O2,$OK2) = gzipwrite($blInf."geneAss","gene cat file [out]",1) if ($reportGeneCat);
	if ($reportEggMapp){
		system "mkdir -p $tmpD" unless (-d $tmpD);
		foreach (my $j=0;$j<@kgdOpts;$j++){
			my $kk = $kgdOpts[$j];
			next if ($kk eq "");
			
			$emFilesFinal{$kk} = "$blInf.emap.$kgdNameShrt[$j].table";
			$emFiles{$kk} = "$tmpD/$DBmode.emap.$kgdNameShrt[$j].table";
			print "$tmpD/$DBmode.emap.$kgdNameShrt[$j].table\n";
			open $OEM{$kk},">$emFiles{$kk}" or die "Can't open $tmpD/$DBmode.emap.$kgdNameShrt[$j].table\n";
	}}
	#die;
	
	my @splSave;
	my @wordv1=(); my @wordv2=();

	while (1){
		my @splNext;
		while (my $line = <$I>){
			chomp $line; 
			my @splX = split (/\t/,$line);
			if (@splX != 12 ){ die "something wrong with blast string!: \n$line\n@splX\n";}
			my $query = $splX [0];
			$query =~ s/\/[12]$// if ($query =~ m/\/[12]/);
			#print $query."\n";
#			$lcnt++;
			
			if ($lastQ eq ""){$lastQ = $query;
			} elsif ($lastQ ne $query){
				$lastQ = $query;
				$stopInMiddle=1;
				@splNext = @splX;
				last;
			}
			
			if ( $splX[0] =~ m/\/2$/){push(@wordv2, [@splX]);
			} else {push(@wordv1, [@splX]);}

		}

		#print "YY ". @wordv2 ." YY ". @wordv1.  "\n";
		#1st combine scores
		#if ($blRes[0]->[0] =~ m/_982/){print "@{$blRes[0]}\n@{$blRes[3]}\n".Dumper(\%wordv1).Dumper(\%wordv2);}
		my $whX = combineBlasts(\@wordv1,\@wordv2);
		#DEBUG:
		#my $whX = \@wordv2;
		#2nd: sort these, find best hit
		#my $arBhit = bestBlHit($whX); #doesn't need this, this is done in the major routine
		#@blRes = @{$whX};#values %{$whX};
		#print "@blRes \n\n";
		#if (exists($wordv1{"GI:BAE78083.1"})){die @blRes."@{$blRes[0]}\n".Dumper(\%wordv1).Dumper(\%wordv2);} #[0]->[0] =~ m/_982/


		#die @blRes."\n";
		#print "Read Assignments..\n";
		for (my $i=0; $i<@aminBLE ; $i++){
			main($whX,$aminBLE[$i],$aminPID[$i],$i,$reportEggMapp);
		}
		undef @wordv1;undef @wordv2; #undef @blRes; 
		if ( $splNext[0] =~ m/\/2$/){push(@wordv2, \@splNext);
		} else {push(@wordv1, \@splNext);}

		
		if ($stopInMiddle==0){last;}
		$lcnt=0;$stopInMiddle=0;
	}
	close $I;
	close $O2 if ($reportGeneCat);
	
	if ($reportEggMapp){foreach my $kk (keys %OEM){	close $OEM{$kk}}}
	
	if ($writeFastaOut){
		foreach my $x (@kgdOpts){
			for (my $i=0; $i<@aminBLE ; $i++){
				close $FAO{$i}{$x};
			}
		}
	}
	if ($writeSumTbls){
		print"writing Tables\n";
		for (my $i=0; $i<@aminBLE ; $i++){
			my $pathXtra = "/CNT_".$aminBLE[$i]."_".$aminPID[$i]."/";
			my $outPath = $inP.$pathXtra;
			writeAllTable($DBmode."parse",$outPath,$aminBLE[$i],$aminPID[$i],$i);#$st1{$i}, $CAThits{$i},$st3{$i});
			#and do eggnog mapping for specific filter parameters
			eggMap_interpret($outPath);
		}
	}
	system "touch $blInf.stone";
	$noMerge /= 2; $mergeDiaHit /= 2;
	print "Merged $mergeDiaHit Diamond hits from DB, $noMerge hits were singletons\n";
	print "all done\n" if ($mode <3);
	
	exit(0);

	
	
	
} elsif ($mode==3) { #only checks for presence of run
	for (my $i=0; $i<@aminBLE ; $i++){
		my $pathXtra = "/CNT_".$aminBLE[$i]."/";
		print "$i  $aminBLE[$i]\n";
		unless(check_files($DBmode."parse",$inP.$pathXtra,$aminBLE[$i])){exit(3);}
	}
} elsif ($mode==4) { #removes run

	for (my $i=0; $i<@aminBLE ; $i++){
		my $pathXtra = "/CNT_".$aminBLE[$i]."/";
		remove_file($DBmode."parse",$inP.$pathXtra,$aminBLE[$i]);
	}
} else {die"unkown run mode!!!\n";}

#print "all done\n";
exit(0);

sub fixCardMeta($){
	my ($hh)=@_; my %dd = %{$hh};
	foreach my $k (keys %dd){
		#print "$dd{$k}->[1]    ";
		$dd{$k}->[1] =~ s/\s/\_/g;
		$dd{$k}->[1] =~ s/-/\./g;
		$dd{$k}->[1] =~ s/[^A-Za-z0-9\.\_\(\)]//g;
		#print "$dd{$k}->[1]\n";
	}
	return (\%dd);
}

sub splitandadd($ $ $){
	my ($hr,$terms,$shareEffect) = @_;
	return if ($terms eq "" || length($terms)==0);
	my @spl = split /,/,$terms;
	my $val = 1;
	$val = 1/scalar(@spl) if ($shareEffect);
	foreach my $term (@spl){
		${$hr}{$term} += $val;
	}
}

sub writeHashTbl($ $ $){
	my ($outF,$hr,$doSortNum) = @_;
	my %hs = %{$hr};
	if (keys %hs == 0){return;}
	open OH,">$outF" or die "Can't open hash out file $outF\n";
	my @kk;
	if ($doSortNum){
		@kk = sort { $hs{$b} <=> $hs{$a} } keys(%hs);
	} else {
		@kk =  keys(%hs);
	}
	foreach my $k (@kk){
		print OH "$k\t$hs{$k}\n";
	}
	close OH;
}

sub eggMap_interpret($){ #higher level annotations with egg nog mapper
	if (!$reportEggMapp){return;}
	my ($outD) = @_;
	print "EggNOG mapper module..\n";
	die "too many evals chosen" if (@aminBLE > 1 ); #fix that only one eval will be printed...
	my @keys1 = keys %emFiles;
	my $tmpEM = "$tmpD/EMcomb$DBmode.table";
	system "rm -f $tmpEM";
	#run only once on combined file..
	my @emFilesDel;
	foreach my $kk (@keys1){
		push (@emFilesDel,$emFiles{$kk});
		system "cat $emFiles{$kk} >> $tmpEM;";
		system "echo 'split_falk_EM_123\t634498.mru_0149\t0\t1000' >> $tmpEM";
		die "$kk entry in emFilesFinal doesn't exist\n" unless (exists($emFilesFinal{$kk}));
		system "rm -f $emFiles{$kk}";
		#last;
	}
	my $oFil = "$tmpD/combEM$DBmode.res"; $oFil =~ s/table//;
	my $cmd = "rm -f $oFil.1 $oFil.2* \n";
	my $KOoutF = "$outD/KOfromEggNogMapper.txt";
	#$print $cmd."\n";
	my %KOperGene;
	if (-e "/g/bork1/huerta/_soft/eggnog-mapper-bigg/emapper.py" && $KOfromNOG){
		$cmd = "/g/bork1/huerta/_soft/eggnog-mapper-bigg/emapper.py --big /g/bork1/huerta/_shared/kegg77_to_egggnog.m8.emapper_table 95 9 -d none -o $oFil.1 --annotate_hits_table $tmpEM";
		print $cmd."\n";
	} else {
		$cmd .= "$emBin -d none --cpu $ncore --no_search --temp_dir $tmpD --no_file_comments --override --no_refine --annotate_hits_table $tmpEM -o $oFil.1\n" ;
	}
	systemW $cmd ;#unless (-e "$oFil.1.emapper.annotations");
	
	#objects to count on higher level
	my %GOs; my %Kmaps; my %COGcats; my %KOs;
	my $ofCnt = 0; my $taxLvlCnt=0;
	#die "@keys1  ".@keys1."\n";
	open O,">$emFilesFinal{$keys1[$ofCnt]}";
	print O "#QueryID\tNOG\tNOGcat\tNOGdescr\tKEGGmap\tGOterms\n";
	open I,"<$oFil.1.emapper.annotations" or die "can't open $oFil.1.emapper.annotations\n";
	while (my $lin = <I>){
		if ($lin =~ m/^#/){next;
		} elsif ($lin =~ m/^split_falk_EM_123\t/){ #split signal to go to next tax lvl
			close O; 
			#write out countups in higher levels
			$emFiles{$taxLvlCnt} =~ m/.emap.([^\.]+).table/;
			$taxLvlCnt++;
			my $taxS  = $1;
			writeHashTbl($outD."eggNogM_${taxS}_GO.txt",\%GOs,1) if (keys(%GOs)>=1);
			writeHashTbl($outD."eggNogM_${taxS}_KEGG_map.txt",\%Kmaps,1) if (keys(%Kmaps)>=1);
			writeHashTbl($outD."eggNogM_${taxS}_COG_cats.txt",\%COGcats,1) if (keys(%COGcats)>=1);
			writeHashTbl($outD."eggNogM_${taxS}_KO.txt",\%KOs,1) if (keys(%KOs)>=1);
			if ($KOfromNOG && keys %KOperGene > 0){
				writeHashTbl($outD."eggNogM_${taxS}_perGene_KO.txt",\%KOperGene,0);
			}
			print "$ofCnt $keys1[$ofCnt] $emFilesFinal{$keys1[$ofCnt]}\n";
			open O,">$emFilesFinal{$keys1[$ofCnt]}" or die "Can't open $emFilesFinal{$keys1[$ofCnt]}\n";
			$ofCnt++;
			next;
		}
		chomp $lin;
		my @spl = split /\t/,$lin;
		my $oStr = "$spl[0]\t";
		my $KOterm = "";
		if ($KOfromNOG && defined $spl[12]){
			my $KOtmp = $spl[12]; $KOtmp =~ s/ko://g;
			$KOperGene{$spl[0]} = $KOtmp;
			$KOterm = $spl[12] ;
		}

		my $Ccat = $spl[10];$Ccat = "" if (!defined $Ccat);
		my $Cdef = $spl[11];$Cdef = "" if (!defined $Cdef);
		my $Kmap = $spl[6]; $Kmap = "" if (!defined $Kmap);
		my $GOterm = $spl[5];$GOterm = "" if (!defined $GOterm);
		if ($spl[8] =~ m/,([^,]+)\@NOG/){
			$oStr .= "\t$1";
			#print $1."\n";
			if ($tabCats==1){
				$Cdef = $COGdef{$1} if (exists($COGdef{$1}));
				$Ccat = $c2CAT{$1} if (exists($c2CAT{$1}));
			} elsif ($tabCats==2){#KEGG
			}
		} else {
			$oStr .= "\t";
		}
		#print "$GOterm , $Kmap , $Ccat\n";
		splitandadd(\%GOs,$GOterm,0);
		splitandadd(\%KOs,$KOterm,0);
		splitandadd(\%Kmaps,$Kmap,1);
		splitandadd(\%COGcats,$Ccat,1);
		$oStr.= "\t$Cdef\t$Ccat\t$Kmap\t$GOterm\n";
		print O "$oStr\n";
	}
	close O; close I;
	system "rm -f $tmpEM $oFil.1.emapper.annotations";

	#write hashes out
	#TODO
	if ($ofCnt < @keys1){die "Not enough taxsplits found in eggNogFile ($ofCnt)\n$oFil.1.annotations\n";}
	#system "rm ".join(" ",@emFilesDel)." $oFil.1*";
}

sub remove_file{
	my ($outF,$outD,$minBLE) = @_;
	my $out = $outD."/".$outF;
	system "rm -f $out*" if (-d $outD);
}
sub check_files{
	my ($outF,$outD,$minBLE) = @_;
	my $out = $outD."/".$outF;
	my @srchFls = (#"$out.$DBmode.GENE2$DBmode", #"$out.GENE2$DBmode"."aCAT"
	"$out.$DBmode.$kgdNameShrt[0].gene.cnts","$out.$DBmode.$kgdNameShrt[0].cat.cnts");
	push (@srchFls,"$out.$kgdNameShrt[0].CATcnts") if ( $DBmode eq "NOG");
	foreach (@srchFls) { unless (-e $_.".gz"){print "Fail $_\n";return 0;}}
	return 1;
}

sub writeAllTable(){
	my ($outF,$outD,$minBLE,$minPID,$pos) = @_;#,$hr1,$hr2,$hr3) = @_;
	#die "X";
	system "mkdir -p $outD" unless (-d $outD);#or die "Failed to create out dir $outD\n";
	my $out = $outD."/".$outF;
	system "rm -f $out*";
	my %COGabundance=%{$COGhits{$pos}}; my %CATabundance=%{$CAThits{$pos}}; my %funHit=%{$GENEhits{$pos}};
	#my %COGabundance=%{$hr1}; my %CATabundance=%{$hr2}; my %funHit=%{$hr3};
	my $cnt = -1;
	foreach my $normMethod (@normMethods){
		$cnt ++;
		foreach my $y (@kgdOpts){
		next if ($kgdName[$y] eq "");
			my @kk; my $O = FileHandle->new;
			if ($tabCats==0){ #MOG CZy etc
				#don't write this for NOGs any longer..
				@kk = sort { $funHit{$normMethod}{$y}{$b} <=> $funHit{$normMethod}{$y}{$a} } keys(%{$funHit{$normMethod}{$y}});
				$O = gzipwrite("$out.$DBmode.$kgdNameShrt[$y].$normMethod.gene.cnts","Dia gene $y Counts");
				#open O,">$out.$DBmode.gene.cnts";
				foreach my $k (@kk){print $O $k."\t".$funHit{$normMethod}{$y}{$k}."\n";}
				close $O;
			}
			if ($tabCats != 2 ){#all but KO & CZy & Patric viru
				@kk = sort { $COGabundance{$normMethod}{$y}{$b} <=> $COGabundance{$normMethod}{$y}{$a} } keys(%{$COGabundance{$normMethod}{$y}});
				$O = gzipwrite("$out.$DBmode.$kgdNameShrt[$y].$normMethod.cat.cnts","Dia cat $y Counts");
				#open O,">$out.$DBmode.cat.cnts"; 
				my $cogCnts=0;
				foreach my $k (@kk){print $O $k."\t".$COGabundance{$normMethod}{$y}{$k}."\n";$cogCnts++;}
				close $O;
				print "Total of $cogCnts categories in $kgdName[$y].\n"if ($cnt==0);
			}
			if ($tabCats){ #NOG || KEGG || CZy || ABRc
				#print CATs 
				@kk = sort { $CATabundance{$normMethod}{$y}{$b} <=> $CATabundance{$normMethod}{$y}{$a} } keys(%{$CATabundance{$normMethod}{$y}});
				#print @kk."\n";
				$O = gzipwrite("$out.$kgdNameShrt[$y].$normMethod.CATcnts","Dia NOG cat $y Counts");
				#open O,">$out.CATcnts";
				my $CATcnt =0;
				foreach my $k (@kk){print $O $k."\t".$CATabundance{$normMethod}{$y}{$k}."\n"; $CATcnt ++;  }#else {print "X";} }
				close $O;
				#print "Total of $CATcnt categories: failed $CATfail of ".($CATfail+$CATexist)." assignments\n";
			}
		}
	}
	#sleep(5);
	#is being counted double due to normethod counts

}

#routine accepts a single read/gene hit to several valid targets
#and selects one hit as "best", annotates this hit & summarizes at higher level
sub main(){
	my ($AR, $minBLE,$minPid,$pos,$reportEggMappNow)=@_;
	my @blRes  = @{$AR};
	if (@blRes == 0){return;}
	#print "@blRes   - ".  @blRes."\n";
	my %COGabundance=%{$COGhits{$pos}}; my %CATabundance=%{$CAThits{$pos}}; my %funHit=%{$GENEhits{$pos}};

#	my %COGabundance=%{$hr1}; my %CATabundance=%{$hr2}; my %funHit=%{$hr3};
	if (!$emode && !$scomode){die "No scoring mode selected!!\n";}
	#print "Scanning $DBmode with an eval of $minBLE\n";
	my $NOGtreat = 0; $NOGtreat = 1 if ($DBmode eq "NOG");
	#my %kgd = %{$kgdHR};
	
	my @tmp = ("");
	@tmp = @{$blRes[0]};
	#print "@tmp PP\n";
	my $qold=$tmp[0];
	die "Subject length not in length DB : $tmp[1]\n" unless (exists ($DBlen{$tmp[1]}));
	my $SbjLen = $DBlen{$tmp[1]};
	my $COGfail=0; my $CATfail=0; my $totalCOG=0; my $CATexist=0;  my $ii=0;
	
	my $eggNOGmap =0;
	if ($minBLE =~ m/^ENM$/){
		$eggNOGmap = 1; #in this mode just print the "best hit" to a file, that egm can work on later..
	}
	#this routine simply counts up number of hits to COGXX
	
	#while (1){ #don't need this, should only go once over this..
	my $fndCat=0;
	
	#my $arBhit = bestBlHit($whX); 

	my $bestSbj="";my $bestID=0; my $bestAlLen=0;my $bestQuery = ""; my $bestIDever=0; my $bestAlLenEver = 0;
	my $COGexists=0; 
	my $bestScore = 0; my $CBMmode = 0; my $bestE=1000; my $bestBitScpre=0;
	my $bestQseq = "";
	#print "XXX ".@blRes."\n";
	for( ;$ii<@blRes;$ii++){
		#print "@{$blRes[$ii]}\n";
		my $Qseq = "";
		my ($Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc) = @{$blRes[$ii]};
		$Qseq = @{$blRes[$ii]}[-1] if ($writeFastaOut);
		#print "$Subject\n";
		#print $Query."  $bitSc\n";
		#sort by eval #changed from bestE -> bestScore
		#die "\n".($cardE{$Subject}**(1/$Card_leniet))."\n$cardE{$Subject}\n";
		#die "Can't find ABRc ID $Subject\n" if ($tabCats==5 && !exists($cardFull{$Subject}) );#|| $cardFull{$Subject}->[2] eq "");
#		if ($Query =~ m/_982/){
#			print "$Query\n$bitSc>=".($cardFull{$Subject}->[3])*($AlLen/$SbjLen)."\n ${$cardFull{$Subject}}[3]  $AlLen /  $SbjLen\n@{$cardFull{$Subject}}\n";
#		}
		next if ($tabCats==5 && !exists($cardFull{$Subject}) );#expected to be missing non relevant 
		if ( ($eval <= $minBLE && $bitSc >= $minScore)
					&& ($AlLen >= $minAlLen)
					&& ($id >= $minPid)
					&& ($quCovFrac == 0 || $AlLen > $SbjLen*$quCovFrac) 
					&& (!$fndCat || $noHardCatCheck || exists $c2CAT{$Subject})
					&& ($tabCats!=5 ||  $bitSc >= ($cardFull{$Subject}->[3]*$AlLen/$SbjLen)  )  #$eval <= ($cardE{$Subject}**(1/$Card_leniet)) ) #Card
					#&& ($tabCats!=5 || $Subject =~ m/^GI:/ )#ABRc specific
					) {
					#now decide if the hit itself is actually better
					#print "BETTER   $c2CAT{$Subject}\n";
#			if (  (($bestID -5)< $id && $bestE*10 > $eval &&
#					( $id >= ($bestID *0.97) && $id >= $bestIDever*0.9   && ( $AlLen >= $bestAlLen * 1.15 ) )  ||  #length is just better (15%+)
#					(   ($id >= $bestID *1.03) && ( $AlLen >= $bestAlLen * 0.9) ) ||#id is just better, while length is not too much off
#					$bitSc > $bestBitScpre*0.8 ) ||  
#					($fndCat==0 && exists $c2CAT{$Subject}) #big one: take in any case
			
			if (  ( ($bestID -5)< $id && $bestE*10 > $eval &&
					( $id >= ($bestID *0.97) && $id >= $bestIDever*0.92   && ( $AlLen >= $bestAlLenEver * 1.15 ) )  ||  #length is just better (15%+)
					(   ($id >= $bestIDever *1.03) && ( $AlLen >= $bestAlLenEver * 0.9) ) ||#id is just better, while length is not too much off
					$bitSc > $bestBitScpre ) ||  
					($fndCat==0 && exists $c2CAT{$Subject}) #big one: take in any case
					){   #convincing score
						#print "ENTER! $id\n";
						$fndCat =1 if (exists $c2CAT{$Subject}); 
						$bestSbj =$Subject; 
						$bestAlLen=$AlLen;$bestE = $eval; $bestQuery = $Query;
						$bestBitScpre=$bitSc;
						$bestID = $id;
						$bestQseq = ">${bestQuery}__$qstart:$qend:$id\n".convertNT2AA($Qseq)."\n" if ($writeFastaOut);
						#die "$bestQseq\n";
						if ($bestID > $bestIDever){$bestIDever=$bestID;}
						if ($bestAlLen > $bestAlLenEver){$bestAlLenEver=$bestAlLen;}
						
			}
		} #else { print "NOT HIT";}
		my @tmp;
	}
#	die;
	#print @blRes ."  $bestSbj $bestQuery $bestID\n";
	#print "$ii ". ($ii+1 >= @blRes)." XX\n";
	#die "$bestQseq\n";
	if ($bestSbj eq ""){return;}
	
	
	#finalize assignment to read
	my $curCat = ""; my $curDef = "";
	my $curCatDef="";
	my $curKgd = 3; #default always 3.. historic
	my @splC;
	if ($ACLdb || $KGMmode){
		@splC = split /:/,$bestSbj ; $splC[0] = $splC[1] if ($ACLdb);
	} elsif ($cazyDB || $tabCats==0 || $tabCats == 6) {
		@splC = split /\|/,$bestSbj;
	} 
	
	##########################################
	#take care of tax assignments / splits
	if ($noTax){
		;
	} elsif ($cazyDB){
		if ($bestSbj =~ m/Bacteria/){$curKgd =0;
		}elsif ($bestSbj =~ m/Archaea/) {$curKgd =1;
		}elsif ($bestSbj =~ m/Eukaryota/){ $curKgd =2 ;
			if (exists($czyTax{$splC[$#splC]})){$curKgd = $czyTax{$splC[$#splC]} ;}
		}else{ #take Mohs long list
			$curKgd = $czyTax{$splC[$#splC]} if (exists($czyTax{$splC[$#splC]}));
		}
		#else {print $bestSbj." ";}
		#get rid of CBM
		if (0 && $bestSbj =~ m/\|CBM\d+/) {$CBMmode = 1; $curKgd = 4;}
	} elsif ( 0 && $ACLdb){ #cats are too primitive, doesn't need to be split up any further
		#67936 >protein:plasmid  25941 >protein:proph  28277 >protein:vir
		if ($bestSbj =~ m/plasmid/){$curKgd =0;
		}elsif ($bestSbj =~ m/proph/){$curKgd =1;
		}elsif ($bestSbj =~ m/vir/){$curKgd =2;
		}
	} elsif ($KGMmode){
		if (exists($KEGGtax{$splC[0]})){
			$curKgd = $KEGGtax{$splC[0]};
		}  #else {print "X";}
	} elsif ($tabCats == 1){ #NOG
		my $taxid = "";
		$bestSbj =~ m/^(\d+)\./; $taxid = $1;
		if (!exists $NOGkingd{$taxid}){print "Can't find $taxid in ref tax\n";}
		$curKgd = $NOGkingd{$taxid} if (!$CBMmode);
	}
	
	#print $curKgd."SSS\n";
	if ($reportEggMappNow){
		print  {$OEM{$curKgd}} "$bestQuery\t$bestSbj\t$bestE\t$bestBitScpre\n" ;
	}
	 my %score;
	

	##############################
	# merged read pair? double the score!
	my $lpCnt=0;
	$score{"cnt"} = 0;
	if ($bestQuery =~ m/\/[12]$/){
		$noMerge++;
		$score{"cnt"}=1;
	} else {
		$mergeDiaHit ++ ; $score{"cnt"}= 2 ;
		#print "HIT";
	}


	
	##############################
	# hierachy assignments etc
	my @curCOGs = ("-");#my $curCOG="-";
	if ($tabCats == 1){
		unless (exists( $g2COG{$bestSbj} )){
			$COGfail++;	return;
		}
		@curCOGs = ($g2COG{$bestSbj});
		$curDef = $COGdef{$curCOGs[0]} if (exists $COGdef{$curCOGs[0]});
	} elsif ($tabCats==7){#AB production
		$bestSbj =~ m/;ID=([0-9\.]+);od=(\S+);ud=(\S+)$/;
		if ($bestID< $1 ){return;}
		@curCOGs = split(/,\s*/,$3);
		#@curCOGs = ($3);
		$curCat = $2;
		my $bestQueryK = $bestQuery; $bestQueryK =~ s/\/\d$//;
		
		if (exists($read2NOGkgd{$bestQueryK})){
			if ( $read2NOGkgd{$bestQueryK} < 2){
				return;
			}
		}
		#die "$bestSbj\n$1\n$2\n$3\n$bestID\n";
	} elsif ($tabCats==5){#ABR CARD
		#$curCOG = $cardName{$bestSbj};
		@curCOGs = ($cardFull{$bestSbj}->[1]);
	}elsif ($tabCats == 2){ #KEGG
		if (exists $c2CAT{$bestSbj} ){
			$curCat = $c2CAT{$bestSbj}; $CATexist++;
		} else {$CATfail++; return;}#print "FAIL: $bestSbj\n\n";
	}elsif ($tabCats==3){  #CAZy
		@curCOGs = split /;/,$splC[0];
	}elsif ($tabCats==6){  #TCDB
		#print "TC6\n@splC\n";
		if ($splC[3] =~ m/(^[^\.]+\.[^\.]+\.[^\.]+)\.?/){
			@curCOGs = ($1);
		} else {
			die "Can't find tag\n";
		}
		if (exists $TCdef{$curCOGs[0]}{CL} ){
			$curCat = $TCdef{$curCOGs[0]}{CL}; $CATexist++;
		} else {$CATfail++; $curCat = $curCOGs[0]; }
		
		$curDef = "";
		if (defined $TCdef{$curCOGs[0]}{desc} ){#CLd1
			$curDef = $TCdef{$curCOGs[0]}{desc}; 
		}
		$curCatDef = "";
		if (defined $TCdef{$curCOGs[0]}{CLd1} ){#CLd1
			$curCatDef = $TCdef{$curCOGs[0]}{CLd1}; 
		}
		#print "$curDef  $curCatDef\n"
		#die $curCOG."\n";
	} elsif($tabCats==8){ #PATRIC virulence
		unless (exists($PTVcat{$bestSbj})){
			die "Can't find PATRICK cat for $bestSbj\n";
		}
		$curCat = $PTVcat{$bestSbj};
		$curDef = $PTVgene_def{$bestSbj};
		@curCOGs = ($PTVgene{$bestSbj});
		#die "$curCOG\n";
	}
	
	if ($writeFastaOut ){
		die "writeFastaOut: key $pos $curKgd not found\n" unless (exists($FAO{$pos}{$curKgd}));
		my @preSeq = split /\n/,$bestQseq;
		$preSeq[0] .= ":$curCat";
		
		print {$FAO{$pos}{$curKgd}} join("\n",@preSeq)."\n";
	}

	#all info parsed, now add up matrices
	my $curCOG =  $curCOGs[0];
	$score{"GLN"} = $bestAlLen / $SbjLen;# if ($normMethod eq "GLN"); #score for this hit, normed by prot length
	my $numCats = scalar @curCOGs;
	foreach my $normMethod (@normMethods){
		my $score2 = $score{$normMethod};
		$lpCnt++;
		$totalCOG += ($numCats*$score{"cnt"}) if ($lpCnt == 1);
		#die "can't find length entry for $bestSbj\n" unless (exists $DBlen{$bestSbj});
		#print $DBlen{$bestSbj}." $normMethod \n";
		if ($tabCats == 1){ #NOG
			#hash{$_} = $valA for qw(a b c);
			#if (!exists($COGabundance{$normMethod}{$curCOG})){$COGabundance{$normMethod}{$curCOG}{$_}=0 foreach (@kgdOpts);}
			$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score2;
			if (exists($c2CAT{$curCOG})){
				$curCat = $c2CAT{$curCOG};$CATexist++;
				#if (!exists($CATabundance{$normMethod}{$curCat})){$CATabundance{$normMethod}{$curCat}{$_}=0 foreach (@kgdOpts);}
				$CATabundance{$normMethod}{$curKgd}{$curCat} += $score2;
			} else {
				$CATfail++;
				#print "Cat unknw for: $curCOG\n" ;
			}
			#print $O "$Query\t$Subject\t$bestID\t$bestE\t$bestAlLen\t$curCOG\t$curCat\t$curDef\n";
			if ($lpCnt==1){
				print $O2 "$bestQuery\t$curCOG\n" if ($reportGeneCat ==1);
				print $O2 "$bestQuery\t$curCOG\t$curCat\t$curDef\n" if ($reportGeneCat ==2);
			}
		
		} elsif ($tabCats == 2 ){ #KEGG
			$CATabundance{$normMethod}{$curKgd}{$curCat} += $score2;
			#print "$bestQuery\t$curCat\n";
			if ($lpCnt==1 && $curCat ne "" ){
				print $O2 "$bestQuery\t$curCat\n" if ($reportGeneCat  );
			}
		} elsif ($tabCats == 8){ # PTV
			$CATabundance{$normMethod}{$curKgd}{$curCat} += $score2;
			$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score2;
			#print "$bestQuery\t$curCat\n";
			if ($lpCnt==1 && $curCat ne "" ){
				print $O2 "$bestQuery\t$curCOG\t$curDef\t$curCat\n" if ($reportGeneCat  );
			}
		} elsif ($tabCats==5){#ABR CARD
			$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score2;
			my @curCats = split(/,\s*/,$cardFull{$bestSbj}->[2]);
			#die "@curCats\n";
			foreach my $cCt (@curCats){
				$CATabundance{$normMethod}{$curKgd}{$cCt}+= $score2;
			}
			if ($lpCnt==1){
				print $O2 "$bestQuery\t$curCOG\n" if ($reportGeneCat ==1);
				print $O2 "$bestQuery\t$curCOG\t".join(",",@curCats)."\n" if ($reportGeneCat ==2);
			}
		} elsif ($tabCats==6){#TCDB
			$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score2;
			$CATabundance{$normMethod}{$curKgd}{$curCat}+= $score2;
			if ($lpCnt==1){
				print $O2 "$bestQuery\t$curCOG\n" if ($reportGeneCat ==1);
				print $O2 "$bestQuery\t$curCOG\t$curDef\t$curCat\t$curCatDef\n" if ($reportGeneCat ==2);
			}
		} elsif ($tabCats==7){#AB production PAB
			#$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score2;
			$CATabundance{$normMethod}{$curKgd}{$curCat}+= $score2;
			#my @curCogs = split(/,\s*/,$curCOG);
			#die "@curCats\n";
			foreach my $cCt (@curCOGs){
				$COGabundance{$normMethod}{$curKgd}{$cCt}+= $score2;
			}

			if ($lpCnt==1){
				print $O2 "$bestQuery\t".join(",",@curCOGs)."\n" if ($reportGeneCat ==1);
				print $O2 "$bestQuery\t".join(",",@curCOGs)."\t$curCat\n" if ($reportGeneCat ==2);
			}
		} elsif ($tabCats==3){  #CAZy
			my @subs2;
			foreach $curCOG (@curCOGs){
				$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score2 ;#/$numCats;
				my @subs = ();
				if (exists ($czySubs{$curCOG})){
					@subs = @{$czySubs{$curCOG}};
					push(@subs2,@subs);
				} else {
					#print "missing cur subject $curCOG\n" ;
				}
				foreach my $cCt (@subs){
					$CATabundance{$normMethod}{$curKgd}{$cCt}+= $score2;#/$numCats;
				}
			}
			#print "@subs\n";
			if ($lpCnt==1){
				print $O2 "$bestQuery\t".join(";",@curCOGs)."\n" if ($reportGeneCat == 1 );
				print $O2 "$bestQuery\t".join(";",@curCOGs)."\t".join(",",@subs2)."\n" if ($reportGeneCat == 2 );
			}
			
		} else { #Moh / MOH2 / CAZY / ACL
			$curCOG = $splC[0];
			$COGabundance{$normMethod}{$curKgd}{$curCOG}+= $score2;
			if ($lpCnt==1){
				print $O2 "$bestQuery\t$curCOG\n" if ($reportGeneCat );
			}
		}
		#if (!exists($funHit{$normMethod}{$bestSbj})){$funHit{$normMethod}{$bestSbj}{$_}=0 foreach (@kgdOpts); } 
		if ($tabCats == 0){ #don't need this for KEGG and NOG
			$funHit{$normMethod}{$curKgd}{$bestSbj}+= $score2;
		}
	}
	
#		$bestSbj="";$bestE=1000; $bestScore=0;$CBMmode=0;$bestBitScpre=0;
#		$bestSbj="";$bestE=10; $bestScore=0;$CBMmode=0;
	$totalCOG /= 2; $CATfail /= 2;#$COGfail /= 2;
	#if (0&& $DBmode eq "NOG"){
	#	print "Total entries:$totalCOG\nCOG not found: $COGfail, CAT not found: $CATfail\n";
	#} 
	#die;
	#return (\%COGabundance, \%CATabundance, \%funHit);
	$COGhits{$pos} = \%COGabundance; $CAThits{$pos} = \%CATabundance;  $GENEhits{$pos} = \%funHit;

}


sub readCOGdef(){
	my ($inF) = @_;
	open I,"<$inF" or die "Can't open $inF\n";
	my %ret;
	while (my $line = <I>){
		chomp $line;
		my ($txx,$COG,$n1,$n2,$cat,$def) = split(/\t/,$line);
		$ret{$COG} = $def;
	}
	close I;
	return \%ret;
}
sub readTCDBdef($){
	my ($inF) = @_;
	#die;
	print "$inF\n";
	open I,"<$inF" or die "Can't open $inF\n";
	my %ret;
	my $newCat=0; my $curClass; my $curClassDesc1;my $curClassDesc2;
	while (my $line = <I>){
		chomp $line;
		unless ($line =~ m/\S/){
			
			$newCat = 1; $curClassDesc1="";$curClassDesc2="";
			next;
		}
		my @spl = split(/ - /,$line);
		if ($newCat){
			$newCat = 0;
			$curClass = $spl[0];
			if (@spl>0){
				$curClassDesc1=$spl[1];
			}
			if (@spl>1){
				$curClassDesc2=$spl[2];
			}
			#die "$curClass  $curClassDesc1  $curClassDesc2\n";
			next;
		}
		my $curGene = $spl[0];
		$ret{$curGene}{CL} = $curClass;
		$ret{$curGene}{CLd1} = $curClassDesc1;
		$ret{$curGene}{CLd2} = $curClassDesc2;
		my $desc = "";
		$desc = $spl[1] if (@spl>0);
		$ret{$curGene}{desc} = $desc;
	}
	close I;
	return \%ret;
}
sub readKeggTax{
	my ($inF) = @_;
	my %ret; my $eukCnt=0; my $proCnt=0;
	open I,"<$inF" or die "Can't open $inF\n";
	while (<I>){
		chomp; my @spl = split /\t/;
		my $taxI = 0; 
		if ($spl[1] eq "euka"){$taxI = 1 ; $eukCnt++} else {$proCnt++;}
		$ret{$spl[0]} = $taxI;
	}
	close I;
	print "Euk in KEGG DB = $eukCnt; Pro = $proCnt\n";
	return \%ret;
}
sub readGene2KO(){
	my ($inF) = @_;
	my %c2cat;
	open I,"<$inF"or die "Can't open $inF\n";
	print "reading gene -> KO file.. " or die "gene2ko file not present\n";
	while (my $line=<I>){
		chomp $line;
		my ($gn,$KO) = split(/\t/,$line);
		$KO =~ s/^ko://;
		$c2cat{$gn} = $KO;
	}
	print "Done\n";
	close I;
	return (\%c2cat);
}


sub readMohTax(){
	my ($inf) = @_;
	my %ret;
	open I,"<$inf" or die "Could not find Moh tax file at $inf\n";
	while (my $line = <I>){
		chomp $line;
		my @spl = split /\t/,$line;
		
		my $taxI =$czyEuk;#euk code
		$taxI = $czyFungi if ($spl[1] =~ m/fungi/i);
		$ret{$spl[0]} = $taxI;
	}
	close I;
	return \%ret;
}
sub createNOGKgd($){
	my ($diaIn) = @_;
	my %rd2Ng;
	my $diaInBest = "$diaIn.Kgd";
	if (-e $diaInBest.".gz"){
		print "Reading pre calc NOG KGD DB\n";
		my $hr = readTabbed3($diaInBest,1);
		%rd2Ng = %{$hr};
		return %rd2Ng;
	}
	#die "Wrong\n";
	#print "$NOGtaxf";
	my %NOGkingd = readNogKingdom($NOGtaxf);
	my ($IE,$OK);
	($IE,$OK) = gzipopen($diaIn,"diamond output file [in]",1); 
	my @splX;
	my $lastQ = ""; my $bEval= 1; my $bSubj = "";
	my $cntKgdHits = 0;
	while ( my $line = <$IE> ) {
		chomp $line; 
		@splX = split (/\t/,$line);
		my $query = $splX[0];
		$query =~ s/\/\d$//;
		if ($lastQ eq ""){$lastQ = $query; 
		} elsif ($lastQ ne $query){
			$bEval = 1;
			#print "$query   $bSubj\n";
			$bSubj =~ m/^(\d+)\./; my $taxid = $1;
			if (!exists $NOGkingd{$taxid}){die "Can't find $taxid in NOG ref tax\n";}
			$rd2Ng{$lastQ} = $NOGkingd{$taxid} ;
			$lastQ = $query; $bSubj = "";
			$cntKgdHits++;
		}
		#print "$bEval > $splX[10]\n";
		if ($bEval > $splX[10]){
			$bEval = $splX[10]; $bSubj = $splX[1];
		}

	}
	close $IE;
	$bSubj =~ m/^(\d+)\./; my $taxid = $1;
	if (!exists $NOGkingd{$taxid}){print "Can't find $taxid in NOG ref tax\n";}
	$rd2Ng{$lastQ} = $NOGkingd{$taxid} ; $cntKgdHits++;

	print "Created new NOG Kgd Hits DB, found $cntKgdHits hits\n";
	my ($OE,$OK2) = gzipwrite($diaInBest,"Best NOG hit file",1);
	foreach my $k (keys %rd2Ng){
		print $OE "$k\t$rd2Ng{$k}\n";
	}
	close $OE;
	#and save this table
	return %rd2Ng;
}
sub readNogKingdom($){
	my ($inT) = @_;
	my %ret;
	open I,"<$inT" or die "Can't open NOG tax file $inT\n";;
	while (my $line = <I>){
		my $cls = 0; #bacteria
		$cls = 1 if ($line =~ m/Archaea/);
		if ($line =~ m/Eukaryota/){
			$cls = 2;
			$cls =3 if ($line =~ m/Fungi/);
		}
		$line =~ m/^(\d+)\s/;
		$ret {$1} = $cls;
	}
	close I;
	return %ret;
}
sub readCzySubs($){#cazy_substrate_info.txt
	my ($inF) = @_;
	open I,"<$inF" or die "Can't open $inF\n";
	my $cnt =0;
	my %ret;
	my $catsInCzu =0 ;
	my @cats = ("PlantCW","Chitin","AlpGlucans","AnimalCarb","BacterialCW","Fructans","FungalCarbs","Dextran");
	while (<I>){
		$cnt++;
		chomp;
		next if ($cnt==1);
		#	Plant Cell Wall Carbohydrates	Chitin	Alpha-glucans	Animal Carbohydrates	Bacterial Cell Wall Carbohydrates	Fructans	Fungal Carbohydrates	Dextran
		my @entr = split /\t/;
		my $idx = shift @entr;
		my $idx2 = $idx; $idx2 =~ s/Unclassified-// ;
		$idx = $idx2 unless (exists($ret{$idx2}));
		#print "x$idx ";
		my $hits=0;
		for( my $x=0;$x< @entr; $x++){
			if ($entr[$x] eq "YES"){
				push (@{$ret{$idx}}, $cats[$x]);
				$hits++;
				$catsInCzu++;
			}
		}
		if ($hits==0){
			$ret{$idx} = [];
		}
	}
	close I;
	#die "$catsInCzu\n";
	return \%ret;
}

sub combineBlasts($ $){
	my ($wh1,$wh2) = @_;
	return $wh1 if (@{$wh2} == 0);
	
	my @Pbl1 = @{$wh1}; my @Pbl2 = @{$wh2};
	#basic dereplication procedure to only take first hit per read
	my %bl1; my %bl2;
	foreach my $li (@Pbl1){
		$bl1{${$li}[1]} = $li unless (exists($bl1{${$li}[1]})); #only take first hit
	}
	foreach my $li (@Pbl2){
		$bl2{${$li}[1]} = $li unless (exists($bl2{${$li}[1]})); #only take first hit
	}
	
	my @ret;
	
	my @allKs = uniq ( keys %bl1, keys %bl2); #
	#print "@allKs\n";
	foreach my $k (@allKs){
		unless (exists ($bl1{$k}) && exists ($bl2{$k}) ){
			next;
			if (exists ($bl1{$k})) {push (@ret,$bl1{$k});#$ret{$k} = $bl1{$k};
			} else { push (@ret,$bl2{$k});}#$ret{$k} = $bl2{$k};}
		}
#		die "$k" unless (exists $bl2{$k});
		my @hit1 = @{$bl1{$k}};
		my @hit2 = @{$bl2{$k}};
		die "Unequal blast string length: \n@hit2 (".@hit2.") \n@hit1(".@hit1.")\n" if (@hit1 != @hit2);
		#pair
		#print "pair";
		#$k =~ s/\/\d$/\/12/;
		my @cur = @hit1;
		$cur[0] =~ s/\/[12]$//;
		
		my @sbss1 = sort{ $a <=> $b }(($hit1[8],$hit1[9])); my @sbss2 = sort{ $a <=> $b }($hit2[8],$hit2[9]); #sort start/stop
		#sort 
		#print "$sbss1[0] > $sbss2[0]\n";
		if ($sbss1[0] > $sbss2[0]){ #sort reads: read1 should start earlier in query than read2
			#print"X $sbss1[0] > $sbss2[0]\n";
			my @tmp = @hit1; @hit1 = @hit2; @hit2 = @tmp;
			@tmp = @sbss1; @sbss1 = @sbss2; @sbss2 = @tmp;
		}
		
		
				#push (@ret,\@hit1); next; #DEBUG

		#query start/stop sorting..
		my @quss1 = sort { $a <=> $b }($hit1[6],$hit1[7]);my @quss2 = sort { $a <=> $b }($hit2[6],$hit2[7]);
		#overlap?
		my $overlap = 0;
		if ($sbss1[1] > $sbss2[0]){ $overlap= ( $sbss1[1] - $sbss2[0]);}# print "OVER $overlap\n";}#die "$cur[3] = $hit1[3] + $hit2[3] - ( $sbss1[1] - $sbss2[0])\n";}
		$cur[3] = ($hit1[3] + $hit2[3]) - $overlap; #ALlength
		if ($cur[3] <0){$cur[3]=$hit1[3];} #fatal.. recover with half useful value
		
		#print "$sbss1[1] > $sbss2[0] $sbss1[0]  $cur[3] $overlap\n";
		$cur[2] = ($hit1[2] + $hit2[2] ) /2.0;#%id
		$cur[11] = ($hit1[11] + $hit2[11]) * (1- $overlap/($hit1[3] + $hit2[3] ) );#bitscore
		$cur[10] = ($hit1[10], $hit2[10])[$hit1[10] > $hit2[10]];  # min($hit1[10] + $hit2[10]); #eval
		#die "@cur\n";
		push (@ret,\@cur);
		#print "$k \n@sbss1 @sbss2\n@hit1\n";
	}
	#my @ret2 = values %ret;
	return \@ret;
}

sub bestBlHit($){
	my ($hr) = @_;
	#my $emode=0; my $scomode = 1;
	my %blasts = %{$hr}; my $bestBit=0; my $bestk=""; my $bestScore=0;my $bestID=0; my $bestLen=0; my $bestIDever=0;
	my $k = "";
	foreach$k(keys %blasts){
		#my ($Query,$Subject,$id,$AlLen,$mistmatches,$gapOpe,$qstart,$qend,$sstart,$send,$eval,$bitSc) = ${$blRes{$k}}[11];
		#print $eval."\n";
		#sort by eval #changed from bestE -> bestScore
		#if ( ( ($emode && $bestE > $eval) || ($scomode && $bestScore< ($id * $AlLen)) ) 
		#			#&& ($eval <= $minBLE || $bitSc >= $minScore)
		#			&& ($quCovFrac == 0 || $AlLen > $DBlen{$Subject}*$quCovFrac) 
		#			&& (($fndCat || $noHardCatCheck) || exists $c2CAT{$Subject}) ) {
		#			#print "Y";
		
		#my $curScore = ${$blasts{$k}}[2] * ${$blasts{$k}}[3];
		#if ($bestScore < $curScore){$curScore = $bestScore;$bestk = $k;	}
		
		#	$bestAlLen=$AlLen;$bestE = $eval;#$bestQuery = $Query;
		
		#just sort by bitscore
		#print "@{$blasts{$k}}\n";
		#if ($bestBit < ${$blasts{$k}}[11]){
		#	$bestBit = ${$blasts{$k}}[11];  $bestk = $k;
		#}
		my $cID = ${$blasts{$k}}[2]; my $cLe = ${$blasts{$k}}[3]; my $cSc=${$blasts{$k}}[11];
		if (($bestID -5)< $cID){
			if ( ( $cID >= ($bestID *0.97) && $cID >= $bestIDever*0.9   && ( $cLe >= $bestLen * 1.15 ) )  ||  #length is just better (15%+)
					(   ($cID >= $bestID *1.03) && ( $cLe >= $bestLen * 0.9) ) ||#id is just better, while length is not too much off
					$cSc > $bestBit*0.8){   #convincing score
				$bestID = $cID;  $bestk = $k; $bestLen = $cLe; $bestBit = $cSc;
				if ($bestID > $bestIDever){$bestIDever=$bestID;}
			}
		}
	}
	if ($bestk eq ""){die "something went wrong with ABR blast:\n$k:$blasts{$k}\n";}
	return $blasts{$bestk};
}

sub help(){
print "Routine to interpret diamond output files to specific datbases (CAZy, KEGG, eggNOG)\n";
print "-i [diamond output]\n";
}





