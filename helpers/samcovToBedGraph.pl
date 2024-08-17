#!/usr/bin/env perl
# A script to convert samtools (>v1.3) depth (with -aa option) output 
# to bedGraph by merging continous positions with the same value.
#https://gist.github.com/taoliu/a44244f98ae2ac0d37282f15167cf69a

#if ( $#ARGV < 0 ) {    print "Need 1 parameter! $0 <samtools depth output (with -aa)>\n";    exit ();}
#$f = $ARGV[ 0 ];
#open( IN,$f);

$_ = <STDIN>;
chomp;
( $chr,$s,$v ) = split;
$e = $s;
$s -= 1;
while( <STDIN> ){
    chomp;
    ( $c_chr,$c_s,$c_v ) = split;
    if( $c_chr eq $chr && $c_v != $v ) {
	print join( "\t",( $chr, $s, $e, $v ) ),"\n";
	$s = $c_s - 1;
	$v = $c_v;
	$e = $c_s }
    elsif( $c_chr eq $chr && $c_v == $v ) {
	$e = $c_s;
    }
    elsif( $c_chr ne $chr ) {
	print join( "\t",( $chr,$s,$e,$v ) ),"\n";
	$chr = $c_chr;
	$s = $c_s-1;
	$v = $c_v;
	$e = $c_s;
    }
}