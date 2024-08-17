#!/usr/bin/env perl
#makes an abundance matrix of NOGs across samples
#./combine_ABR.pl /g/scb/bork/hildebra/SNP/GNMass2_singl/
#./combine_ABR.pl /g/scb/bork/hildebra/Tamoc/FinSoil/
use warnings;
use strict;
use Mods::GenoMetaAss qw(gzipopen readMap);
sub writeMat;

my @DBs = ("ABR"); # ("KGB","KGE","CZy","NOG","ABR","MOH");#("NOG","CZy","MOH");
my $rareBin = "/g/bork3/home/hildebra/dev/C++/rare/./rare";

my $inD = $ARGV[0];
$inD .= "/" unless ($inD =~ m/\/$/);
my $outD1 = $inD."pseudoGC/FUNCT/ABR_k/";
system "mkdir -p $outD1" unless (-d $outD1);
my ($hrm,$hr2) = readMap($inD."LOGandSUB/inmap.txt",0,{},{},0);
my %map = %{$hrm};
my @samples = @{$map{opt}{smpl_order}};
my %cat;


#look in each subdir specific to DB for valid *cat file to get the TAXs available
foreach my $DB (@DBs){
	foreach my $NORM ("cnt"){#,"GLN"){
		foreach my $smpl(@samples){
			print "$smpl  XX\n";
			my $testD = $map{$smpl}{wrdir}."/diamond/ABR/";
			my $ABRfile = $testD."ABR.cats.txt";
			next unless (-e $ABRfile);
			open I,"<$ABRfile" or die "couldn't find $ABRfile\n";
			while (my $line = <I>){
#				print $line;
				my $val = 2;#add 1 or 2 read counts?
				chomp $line; 
				my @spl = split /\t/,$line;
				$val = 1 if ($spl[0] =~ m/\/\d$/);
				#print "$spl[0] $val\n";
				for (my $sk=0;$sk<3;$sk++){
					next if (@spl <= 3+$sk);
					my @spl2 = split /,/,$spl[3+$sk];
#					print $spl[3+$sk]." ";
					foreach my $sps (@spl2){
						$cat{$sk}{$sps}{$smpl} += $val;
					}#sps
				}
			}#read lines
			close I;
		}#samples
	}#NORM
}#DB

#write matrices out
for (my $sk=0;$sk<3;$sk++){
	open O,">$outD1/ABRcat.$sk.txt";
	print O "CAT$sk\t".join ("\t",@samples) . "\n";
	foreach my $sps (keys %{$cat{$sk}}){
		print O $sps;
		foreach my $smpl(@samples){
			if (exists $cat{$sk}{$sps}{$smpl}){print O "\t".$cat{$sk}{$sps}{$smpl};
			} else {print O "\t0";}
		}
		print O "\n";
	}
	close O;
}

print "All done\n";
exit(0);
