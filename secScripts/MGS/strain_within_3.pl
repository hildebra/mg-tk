#!/usr/bin/perl
#uses within species (MGS) trees to select subsets of samples that closely related These will be aligned at much higher resolution
#actually original idea is not relevant, what I need is families / individuals aligned at high res, with ~100 samples max
use warnings;
use strict;

use Mods::GenoMetaAss qw( readClstrRev systemW readMapS readFasta);
use Mods::Subm qw(qsubSystem emptyQsubOpt);
use Mods::IO_Tamoc_progs qw(getProgPaths );
use Mods::geneCat qw(readGene2tax createGene2MGS);
my $treeSubGrpsR = getProgPaths("treeSubGrpsR");
