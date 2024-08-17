#!/usr/bin/env perl

use warnings; use strict;

#http://pfam.xfam.org/family/Transpeptidase#tabview=tab6
my $hmmBin3 = "/g/bork3/home/hildebra/bin/hmmer-3.1b1-linux-intel-x86_64/binaries/hmmsearch";

my $hmmModel = "/g/scb/bork/hildebra/DB/HMMs/Transpeptidase.hmm";
my $faa = "/g/bork5/hildebra/results/TEC2/6666666.162524.faa";
print "$hmmBin3 -Z 11927849 -E 1000 --cpu 4 $hmmModel $faa >/g/bork5/hildebra/results/TEC2/TransPep.hmm.res\n";

print "samtools faidx $faa 'fig|6666666.162524.peg.1641' 'fig|6666666.162524.peg.2390' >/g/bork5/hildebra/results/TEC2/TransPep.hmm.res.faa\n";