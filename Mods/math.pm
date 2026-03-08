package Mods::math;
use warnings;
use Cwd 'abs_path';
use strict;
#use List::MoreUtils 'first_index'; 
use Mods::IO_Tamoc_progs qw(getProgPaths);

use Exporter qw(import);
our @EXPORT_OK = qw(avgArray medianArray quantileArrayR nonZero meanArray roundAr round roundF quantileArray); 


#N non-zero values in array
sub nonZero($){
	my ($ar) = @_;
	my @a = @{$ar};
	my $nc =0 ;
	foreach my $x (@a){ $nc ++ if ($x > 0);}
	return $nc;
}

#simply calcs the average array
sub avgArray{
	my ($AR) = @_;
	my @AoA = @{$AR};
	my $n = scalar @AoA;
	if ($n == 1){return \@{$AoA[0]};}
	
	my @ret = @{$AoA[0]};
	for (my $i=0;$i<@ret;$i++){$ret[$i] /= $n;} #norm first entry
	for (my $j=1;$j<$n;$j++){
		for (my $i=0;$i<@ret;$i++){
			$ret[$i] += ${$AoA[$j]}[$i] / $n;
		}
	}
	return \@ret;
}

sub roundAr{
	my ($AR,$ndp) = @_;
	my $factor = 10**$ndp;
	my @A = @{$AR};
	
	for (my $i=0;$i<@A;$i++){
		$A[$i] = int($A[$i] * $factor+0.5)/$factor;
	}
	return \@A;
}

sub round{
	my ($val,$ndp) = @_;
	my $factor = 10**$ndp;
	return (int($val * $factor+0.5)/$factor);
}
sub roundF{
	my ($val,$factor) = @_;
	return (int($val * $factor+0.5)/$factor);
}
sub meanArray{
	return if (@_ == 0);
	my @AR = @{$_[0]};
	my $n = scalar @AR; my $sa=0;
	foreach my $x (@AR){
		$sa += $x;
	}
	return ($sa/$n);
}
sub medianArray
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {        return $vals[int($len/2)];
    }    else #even
    {        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub quantileArray
{
	my $frac = shift;
	return 1 if (scalar(@_) == 0);
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    return $vals[int($len*$frac)];
}

sub quantileArrayR {
    # Usage:
    #   my $q50  = quantileArray(\@arr, 0.5);
    #   my ($q10, $q90) = quantileArray(\@arr, 0.1, 0.9);

    my ($aref, @thresholds) = @_;

    die "quantileArray: need an array ref and at least one threshold\n"
        unless ref($aref) eq 'ARRAY' && @thresholds;
    die "quantileArray: array is empty\n"
        unless @$aref;

    for my $t (@thresholds) {
        die "quantileArray: threshold $t out of range [0,1]\n"
            unless $t >= 0 && $t <= 1;
    }

    my @sorted = sort { $a <=> $b } @$aref;
    my $n      = scalar @sorted;

    my @results;
    for my $t (@thresholds) {
        my $pos = $t * ($n - 1);      # floating point position in sorted array
        my $lo  = int($pos);          # lower index
        my $hi  = $lo + 1;            # upper index
        my $frac = $pos - $lo;        # interpolation fraction

        my $quantile;
        if ($hi >= $n) {
            # $t == 1.0 edge case: clamp to last element
            $quantile = $sorted[-1];
        } else {
            # linear interpolation between neighbours
            $quantile = $sorted[$lo] * (1 - $frac) + $sorted[$hi] * $frac;
        }
        push @results, $quantile;
    }

    return wantarray ? @results : $results[0];
}