#! /usr/bin/perl -w
#perl site_taxon_table.pl Family * > family.ITS.tab ### usage

use strict;
use Mods::GenoMetaAss qw(gzipopen);

my $tax_level_a= lc shift @ARGV;
my $outF = shift @ARGV;
my @tlvls = split(/,/,$tax_level_a);
my $inDir = shift @ARGV;
my %sites;
my %taxa;
my $numLvl = scalar @tlvls;

opendir(DIR, $inDir) || die "can't opendir $inDir: $!";
my @files = grep { /\.hiera\.txt\.gz/ && -f "$inDir/$_" } readdir(DIR);
closedir DIR;

print "Detected ".@files." input files in dir $inDir\n";
exit(0) if (@files ==0);
my %column ;

foreach my $file (@files) {
		#open my $FHANDLE, "< $file" or die "Could not open $file: $!\n";
		my ($FHANDLE,$readinOk) = gzipopen("$inDir/$file","tax infile");
		my $tag = $file;
		$tag =~ s/.*\///;
		$tag =~ s/\.hiera\.txt$//;
         my $cset = 0;
        while (my $row=<$FHANDLE>) {
			chomp $row;
			
			my @temp = split /\t/, $row;
			if ($cset == 0) {
				foreach my $l (@tlvls){
					for (my $i=0; $i<scalar @temp; ++$i) {
							if (lc $temp[$i] eq $l) { $column{$l} = ($i-1); last; } #new LCA: get rid of first entry
							if ($i+1 == scalar @temp) { die "Could not find given taxon level $l in $file\n"; }
					}
					#$sites{$l}{$tag} = {};
				}
				$cset=1;
			} else {
				shift @temp; #new LCA: get rid of first entry
				#foreach my $l (@tlvls){
				#rm Opisthokonta from PR2 DB.. annoying
				if ($temp[1] eq "Opisthokonta"){
					splice @temp, 1, 1;
					splice @temp, 4, 0, "?";
				}
				for (my $ii=0;$ii< $numLvl; $ii++){
					#die "\n".$column+1 ."\n";
					my $k = "";
					#if (@temp <= $column{$l}){
					if (@temp <= $ii){
						$k = join (';',@temp) . join(";", "?" x ($ii - $#temp));
					} else {
						$k = join (';',@temp[0 .. $ii]);
					}
					#print "$column{$l} @temp\n$k\n$file\n";
					#die $k."\n";
					++$sites{$ii}{$tag}{$k};
					++$taxa{$ii}{$k};
				}
			}
		}
		close $FHANDLE;
}
print "Read input files..\n";
foreach my $l (@tlvls){
	my $ii = $column{$l};
	my @taxa_keys = sort {$taxa{$ii}{$b} <=> $taxa{$ii}{$a}} keys %{$taxa{$ii}};
	my @sites_keys = sort keys %{$sites{$ii}};
	my %locSites = %{$sites{$ii}};
	open O,">$outF.$l.txt"; print O "$l";
	foreach my $site (@sites_keys) { print O "\t$site"; }
	print O "\n";
	foreach my $key (@taxa_keys) {
		print O $key;
		foreach my $site (@sites_keys) { 
			if (exists($locSites{$site}{$key})) { print O "\t$locSites{$site}{$key}"; }
			else { print O "\t0"; }
		}
		print O "\n";
	}
	close O;
	system "gzip $outF.$l.txt";
}