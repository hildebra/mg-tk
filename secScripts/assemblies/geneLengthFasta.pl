#!/usr/bin/env perl
use Mods::GenoMetaAss qw(gzipwrite gzipopen readFasta);

my ($inF,$outF) = @ARGV;

my ($I,$OK) = gzipopen($inF,"gene fasta",1);
die "Can't open $inF" unless ($OK);
open O,">$outF" or die "Can;t open output DB length file $outF\n";
my $seqL=0;my $hd="";
while (my $line = <$I>){
	chomp $line;
	#print $line."\n";;
	if ($line =~ m/^>(\S+)/){
		#next if ($seqL==0);
		
		print O "$hd\t".$seqL."\n" unless ($hd eq "");
		$seqL = 0; $hd = $1;
		#print $hd."\n";
		next;
	}
	$seqL += length($line);
}
print O "$hd\t".$seqL."\n" ;


#my $fr = readFasta($inF,1);
#my %FA = %{$fr};
#foreach my $k (keys %FA){
#	print O "$k\t".length($FA{$k})."\n";
#}

close O; close $I;

print "Done\n";