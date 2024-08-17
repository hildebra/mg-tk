#use strict;
#use warnings;
use Mods::TamocFunc qw(sortgzblast);


$ardbfile = "/g/bork1/forslund/tara_resistome/ardb.tabs.parsed";
$mapfile = "/g/bork1/forslund/tara_resistome/ardb_and_reforg_mapping";
$besthitfile = "/g/bork1/forslund/tara_resistome/ardb_vs_reforg9f.overlap90shortest_famthres_or_symbol.sorted.besthit";
$outputfile = $ARGV [1]; # "/g/bork1/forslund/tara_resistome/test.gz";
$outputfilecats = $ARGV [2]; # "/g/bork1/forslund/tara_resistome/testCats.gz";
$inputfile = $ARGV [0]; # "/g/scb/bork/hildebra/Tamoc/FinSoil/Sample_AV110_4/diamond/dia.ABR.blast.gz";

$inputfile = sortgzblast($inputfile);
 
 
$outputfilecats =~ m/^(.*\/)[^\/]+$/;
my $outD = $1;
system "mkdir -p $outD" unless (-d $outD);
 
#reads ardb tabs
%ssym = (); #cat 1
%scat = (); #cat 2
%sthres = (); #cat 3

#read CAT DB
open (FH, $ardbfile);
while (<FH>) {

    $aLine = $_;
    chomp ($aLine);

    @words = split (/\t/, $aLine);

	if (uc ($words [3]) eq "BACA") { $words [7] = "80"; } # assume/fix typo in ardb file
    $ssym {$words [0]}{uc ($words [2])} = "";
    $scat {$words [0]}{uc ($words [3])} = "";
    $sthres {$words [0]}{uc ($words [7])} = "";
}
close (FH);

open (FH, $mapfile);
%sym2drug = (); #specific drug resistance
while (<FH>) {
    $aLine = $_;
    chomp ($aLine);
    @words = split (/\t/, $aLine);
    $sym2drug {$words [1]} = $words [3];
}
close (FH);


#read in ID cutoffs
open (FH, $besthitfile); # besthit file
while (<FH>) {
	$aLine = $_;
	chomp ($aLine);
	@words = split (/\t/, $aLine);
	foreach $asym (keys % {$ssym {$words [1]}}) {
		$ssym {$words [0]}{$asym} = "";
	}
	foreach $acat (keys % {$scat {$words [1]}}) {
		$scat {$words [0]}{$acat} = "";
	}
	foreach $athres (keys % {$sthres {$words [1]}}) {
		$sthres {$words [0]}{$athres} = "";
	}
}

close (FH);

open (FH2, " > $outputfile");
open (FH3, " > $outputfilecats");
open (FH, "zcat $inputfile |");

while (<FH>) {
	$aLine = $_;
	chomp ($aLine);

	@words = split (/\t/, $aLine);
	$query = $words [0];
	$subject = $words [1];
#	$qlen = $words [3];
#	$slen = $words [4];
	$id = $words [2];
#	$bitscore = $words [5];
#	$all = $words [6];
	$evalue = $words [10];

##	$overlap = $all / $llen;
	@thress = sort {$a <=> $b} keys % {$sthres {$subject}};
	$thres = $thress [0];
	if ($evalue < 0.00001 && $id >= $thres) {
		%ldrugs = ();
		foreach $sym (keys % {$ssym {$subject}}) {
			foreach $drug (split (/\,/, $sym2drug {$sym})) {
				$ldrugs {$drug} = "";
			}
		}

		print FH2 "$aLine\n";
		print FH3 "$query\t$subject\t$id\t".join (",", keys % {$ssym {$subject}})."\t".join (",", keys % {$scat {$subject}})."\t".join (",", keys % ldrugs)."\n";
	}
}

close (FH);
close (FH2);
close (FH3);


system "touch $inputfile.stone";
