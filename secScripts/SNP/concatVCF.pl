#!/usr/bin/env perl
#takes several vcf and creates a new vcf with all SNPs merged
#usage: ./concatVCF.pl vcf1 vcf2 ..
use warnings;
use strict;
use Mods::GenoMetaAss qw(readMap );
use Mods::Subm qw(qsubSystem emptyQsubOpt);

my %SNPs;
my $HEAD1=""; my $HEAD2=""; my $HEADref=""; my $HEADsrc="##source="; my $HEADset = 0; 
#some stats
my %multHit; my %uniqHit;

foreach my $ivcf (@ARGV){
my $sourceFnd=0;
	open I,"<$ivcf" or print STDERR "cant' open $ivcf\n";
	while (my $line = <I>){
		chomp $line;
		if ($line =~ m/^#/){
			if ($line =~ m/##source=(.*)$/){
				$HEADsrc .= "$1;";
				$sourceFnd=1;
			} elsif ($HEADset==0){
				if ($sourceFnd==0){$HEAD1.=$line."\n";} else {$HEAD2.=$line."\n";}
			}
			next;
		} elsif ($HEADset==0){$HEADset=1;}
	
		my @spl = split /\t/,$line;
		if (exists($SNPs{$spl[0]}{$spl[1]})){
			#die "$line\n$SNPs{$spl[0]}{$spl[1]}\n";
			$multHit{$spl[0]}++;
		} else {
			$uniqHit{$spl[0]}++;
		}
		$SNPs{$spl[0]}{$spl[1]} = $line;
	}
	close I;
}


print $HEAD1.$HEADsrc."\n".$HEAD2;
foreach my $k (sort (keys %SNPs)){
	foreach my $k2 (sort {$a <=> $b} (keys %{$SNPs{$k}})){
		print $SNPs{$k}{$k2}."\n";
	}
}

#print STDERR "Mult SNP=$multHit; UNIQ SNP=$uniqHit;\n";
print STDERR "MULT / UNIQ SNPs = "; foreach my $k (keys(%multHit)){print STDERR "$k:$multHit{$k}/$uniqHit{$k} ; ";}print STDERR "\n";



