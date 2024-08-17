#!/usr/bin/perl
use strict;
use warnings;

my $inD = $ARGV[0];

opendir(DIR, $inD);
my @files = grep(/\.fna\.depStat$/,readdir(DIR));
closedir(DIR);

my %mat;my@smpls;my %mAlt; my %mRef; my %mFreq;
foreach my $f (@files){
	my $read=0; my $r1=0;
	$f =~ m/.*_(.*)\.cons\.fna\.depStat/;
	my $smpl = $1;
	push @smpls,$smpl;
	open I,"<$inD/$f" or die "Can;t find file $inD/$f\n";
	while (<I>){
		chomp;
		if (m/altFreq	0	0/){$r1=1;next;}
		my @spl = split /\t/;
		if ($r1){
			if (m/Alternate allele freqs:/){$read=1 ; $r1 = 0;next;}
			next if (@spl ==0);
			my $cd = $spl[0];
			if (@spl < 3){die "@spl\n$inD/$f\n";}
			$mAlt{$cd}{$smpl} = $spl[1];
			$mRef{$cd}{$smpl} = $spl[3];
			$mFreq{$cd}{$smpl} = $spl[2];
		}
		next unless ($read);
		
		$mat{$spl[0]}{$smpl} = $spl[1];
	}
	close I;
}
#die;

open O,">$inD/Stat_freq.txt";
foreach my $s (@smpls){print O "\t$s";}
print O "\n";
foreach my $k1 (sort {$a <=> $b} keys %mat){
	print O $k1;
	foreach my $s (@smpls){
		if (exists($mat{$k1}{$s})){	print O "\t$mat{$k1}{$s}";
		} else {print O "\t0";	}
	}print O "\n";
}
close O;
open O,">$inD/Stat_Ref_cnt.txt" or die "Can't open $inD/Stat_Ref_cnt.xtxt\n";
foreach my $s (@smpls){print O "\t$s";}
print O "\n";
foreach my $k1 (sort {$a <=> $b} keys %mRef){
	print O $k1;
	foreach my $s (@smpls){
		if (exists($mRef{$k1}{$s})){	print O "\t$mRef{$k1}{$s}";
		} else {print O "\t0";	}
	}print O "\n";
}
close O;
#die "$inD/Stat_Ref_cnt.txt";
open O,">$inD/Stat_AltFreq.txt";
foreach my $s (@smpls){print O "\t$s";}
print O "\n";
foreach my $k1 (sort {$a <=> $b} keys %mFreq){
	print O $k1;
	foreach my $s (@smpls){
		if (exists($mFreq{$k1}{$s})){	print O "\t$mFreq{$k1}{$s}";
		} else {print O "\t0";	}
	}print O "\n";
}
close O;
open O,">$inD/Stat_Alt_cnt.txt";
foreach my $s (@smpls){print O "\t$s";}
print O "\n";
foreach my $k1 (sort {$a <=> $b} keys %mAlt){
	print O $k1;
	foreach my $s (@smpls){
		if (exists($mAlt{$k1}{$s})){	print O "\t$mAlt{$k1}{$s}";
		} else {print O "\t0";	}
	}print O "\n";
}
close O;

print "$inD/FreqStat.txt\n";
