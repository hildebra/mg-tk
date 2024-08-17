#!/usr/bin/perl -w

#
# (C) 2003 by
#     Markus Goeker (markus.goeker@uni-tuebingen.de) and
#     Clemens Oertel (oertel@uni-tuebingen.de),
#     University of Tuebingen, Germany
#
# This program is distributed under the terms of the Gnu Public License V2.
# For further information, see http://www.gnu.org/licenses/gpl.html
#

##############################################################################
#
# GLOBAL VARIABLES
#
#############################################################################

my $DEBUG = 1;
my $OPT_ONE = 0;
my $OPT_USE_WEIGHTS = 0;

my @AVERAGES = ();
my @STDDEVS = ();
my @VARIANCES = ();
my @VARIANCESINV = ();
my @TRANSLATIONS = ();
my @TREES = ();

my $CURRENT_INLINE = '';
my $EOIF = 0;
my $SECTION = '';
my @SECTION_ADDONS = ();

my $IN = 0;
my $OUT = 0;

my $FILETYPE;
my $GETLINE;
my $OPENFILE;
my $READTAXA;

##############################################################################
#
# FILE IO
#
#############################################################################

sub det_filetype {
    my ($if) = @_;
    
    my $sep = '';
    
    open(FH, "<$if") || die "Cannot open $if, aborting ...\n";
    
    do { $sep = getc(FH); } until ($sep eq "\012" || $sep eq "\015");

    if (($sep eq "\015") && (getc(FH) eq "\012")) {
	$sep .= "\012";
    }
    
    $/ = $sep;
    seek(FH, 0, 0);

    while ((! $FILETYPE) && (my $l = <FH>)) {
	if ($l =~ /^#NEXUS/i) {
	    $FILETYPE = "nexus";
	    
	    $OPENFILE = \&iox_openfile;
	    $GETLINE = \&iox_getline;
	    $READTAXA = \&iox_readtaxa;
	} elsif ($l =~ /^\s*\d+\s*$/) {
	    $FILETYPE = "phylip";
	    
	    $OPENFILE = \&iop_openfile;
	    $GETLINE = \&iop_getline;
	    $READTAXA = \&iop_readtaxa;
	}
    }
    
    die "Could not determine type of $if, aborting ...\n" unless $FILETYPE;
    
    close(FH);
}

sub iox_openfile {
    my ($if) = shift;
    
    open(FH, "<$if") || die "Cannot open $if, aborting ...\n";
    
    $EOIF = 0;
    $CURRENT_LINE = '';
    
    my $header = <FH>;
    die "$if is not a nexus-file, aborting ...\n" unless $header =~ /^#NEXUS/i;

    return *FH;
}

sub iox_getline {
    my $line = {};
    my $comment = 0;
    
    while ($CURRENT_INLINE !~ /;/) {
	if (defined($IN) && (my $l = <$IN>)) {
	    chop $l;
	    
	    next if $l =~ /^\s*$/;
	    next if $l =~ /^#/;

	    # End of multiline comment
	    if ($l =~ /^([^\[]*)\](.*)/) {
		$line->{'addons'}[$#{$line->{'addons'}}] .= ' ' . $1;
		$l = $2;
		
		$comment = 0;
	    }
	    
	    # Within multiline comment
	    if($comment) {
		$line->{'addons'}[$#{$line->{'addons'}}] .= ' ' . $l;
		next;
	    }
	    
	    # Begin of multiline comment
	    if ($l =~ /(.*)\[([^\]]*)$/) {
		$l = $1;
		push(@{$line->{'addons'}}, $2);
		
		$comment = 1;
	    }
	    
	    # One-line/partial-line comments
	    while ($l =~ /\[([^\]]*)\]/) {
		push(@{$line->{'addons'}}, $1);
		$l =~ s/\[[^\]]*\]//;
	    }
	    
	    $CURRENT_INLINE .= ' ' . $l;
	} else {
	    $EOIF = 1;
	    last;
	}
    }
    
    $CURRENT_INLINE =~ s/\s*;\s*/;/g;
    $CURRENT_INLINE =~ s/^\s*//g;
    $CURRENT_INLINE =~ s/\s*$//g;
    $CURRENT_INLINE =~ s/\s+/ /g;
    
    if ($CURRENT_INLINE =~ /^([^;]+);/) {
	$line->{'data'} = $1;
	$CURRENT_INLINE =~ s/^([^;]+);//;
    } else {
	$line->{'data'} = $CURRENT_INLINE;
	$CURRENT_INLINE = '';
    }
    
    if ($line->{'data'} =~ /^begin\s+([^\s]+)$/i) {
	@SECTION_ADDONS = @{$line->{'addons'}} if defined($line->{'addons'});
	$SECTION = lc($1);
    } elsif ($line->{'data'} =~ /^end(block)?$/i) {
	@SECTION_ADDONS = ();
	$SECTION = '';
    } else {
	push(@{$line->{'addons'}}, @SECTION_ADDONS);
    }

    return ($line->{'data'} ? $line : 0);
}

sub iox_readtaxa {

}

sub iop_openfile {
    my ($if) = shift;
    
    open(FH, "<$if") || die "Cannot open $if, aborting ...\n";

    $EOIF = 0;
    $CURRENT_LINE = '';
    
    return *FH;
}

sub iop_getline {
    my $data;
    
    while ((! $CURRENT_LINE) && (! $EOIF)) {
	defined($IN) and $CURRENT_LINE = <$IN> or $EOIF = 1;
    }
    
    while ((! $data) && (! $EOIF)) {
	my $squared = 1;
	my $row = -1;
	my $cnt = 0;
	my $t = 1;
	
	if ($CURRENT_LINE =~ /^\s*(\d+)\s*$/) {
	    my $n = $1;
	    
	    defined($IN) and $CURRENT_LINE = <$IN> or $EOIF = 1;
	    $data = 'matrix';
	    
	    last if $EOIF;

	    do {
		if ($CURRENT_LINE =~ /^([^\s]+)((\s+[e\d\.\-]+)+)$/) {
		    $TRANSLATIONS[$t] = $1; $t++;

		    my $l = $2;
		    $cnt = 0; $row++;
		    
		    if ($squared) {
			while (($l =~ /^\s*([e\d\.\-]+)(.*)$/) && ($cnt < $row)) {
			    $data .= ' ' . $1;
			    $l = $2;
			    $cnt++;
			}

		#	$row++ if $cnt && $cnt == $row;
		    } else {
			$data .= ' ' . $l;
		    }
		} elsif ($CURRENT_LINE =~ /^((\s+[e\d\.\-]+)+)$/) {
		    my $l = $CURRENT_LINE;
		    
		    if ($squared) {
			if ($cnt < $row) {
			    while (($l =~ /^\s*([e\d\.\-]+)(.*)$/) && ($cnt < $row)) {
				$data .= ' ' . $1;
				$l = $2;
				$cnt++;
			    }
			    
		#	    $row++ if $cnt && $cnt == $row;
			}
		    } else {
			$data .= ' ' . $l;
		    }
		} elsif ($CURRENT_LINE =~ /^([^\s]+)$/) {
		    $TRANSLATIONS[$t] = $1; $t++;
		    
		    $cnt = 0; $squared = 0;
		} else {
		    die "Format error in input file ($l), aborting ...\n";
		}

		defined($IN) and $CURRENT_LINE = <$IN> or $EOIF = 1;
		print "LINE: $data\n";
	    } until ($EOIF || ($CURRENT_LINE =~ /^\s*\d+\s*$/));
	} else {
	    die "Format error in input file ($l), aborting ...\n";
	}
    }
    
    return unless $data;
    
    $data =~ s/^\s*//g;
    $data =~ s/\s*$//g;
    $data =~ s/\s+/ /g;
    
    $SECTION = "distances";

    return { 'data' => $data };
}

sub iop_readtaxa {

}

sub iop_phylip2nexus {
    my ($if) = @_;
    
    $IN = &iop_openfile($if);

    $if =~ /^(.*)\.[^\.]+$/;
    open($OUT, ">$1.nex");

    print $OUT "#NEXUS\n";
    print $OUT "[! Distance matrices from phylip-formatted distance file " . $if . "]\n";
    
    my $no = 0;
    my $line = &iop_getline($IN);
    &print_taxa();
    
    do {
	if ($SECTION eq "distances") {
	    if ($line =~ /matrix\s+(([e\d\.\-]+\s*)+)/) {
		my @values = split(/\s+/, $1);
		$no++;
		
	        print $OUT "[! PHYLIP_$no]\n";
	        print $OUT "begin distances;\n";
	        print $OUT "  format triangle=lower nodiagonal nolabels;\n";
	        print $OUT "  matrix\n";
	        &print_ALTmatrix(\@values);
	        print $OUT "  ;\n";
	        print $OUT "end;\n";
	    }
	}
    } while ($line = &iop_getline($IN));
    
    close($OUT); close($IN);
}

##############################################################################
#
# TREES --> DISTANCE MATRICES
#
#############################################################################

sub process_tree_file {
    my ($if, $of, $from, $to) = @_;
    
    my $i = 0;
    
    print "Processing trees ...\n" if $DEBUG;
    
    $IN = &$OPENFILE($if);
    open($OUT, ">$of");
    
    print $OUT "#NEXUS\n";
    print $OUT "[! Patristic distance matrices from treefile " . $if . "]\n";
    
    while (my $line = &$GETLINE($IN)) {
	if ($SECTION eq "trees") {
	    if ($#TRANSLATIONS == -1 && $line->{'data'} =~ /^translate/i) {
		while ($line->{'data'} =~ /translate\s+(\d+)\s+([^,]+)/i) {
		    $TRANSLATIONS[$1] = $2;
		    $line->{'data'} =~ s/translate\s+\d+\s+[^,]+(,\s*)?/translate /i;
		}
		
		&print_taxa();
		next;
	    }

	    if ($line->{'data'} =~ /^u?tree(?:\s+|\s*\*\s*)([^\s]+)\s*=\s*([^;]+)$/i) {
		die "No translations found in $if, aborting ...\n" if $#TRANSLATIONS == -1;
		print "  Tree $1\n" if $DEBUG;

		$i++;
		&process_tree($1, $2, $line->{'addons'}) if $i >= $from && ($to == 0 || $i <= $to);
	    }
	}
    }
    
    close($OUT); close($IN);
    
    print "done.\n" if $DEBUG;
}

sub process_tree {
    my ($name, $tree, $aref) = @_;
    
    my $vid = $#TRANSLATIONS;
    my $tref = { 'name' => $name, 'addons' => $aref };
    
    while ($tree =~ /\(([^\(\)]+)\)/) {
	my $twig = $1;
	
	$vid++;
	$tree =~ s/\(([^\(\)]+)\)/$vid/;
	
	if ($twig =~ /\d:\d/) {
	    while ($twig =~ /^(\d+):([e\d\.\-]+)/) {
		my $id = $1;
		my $dist = ($OPT_ONE ? 1 : $2);
		
		$twig =~ s/^(\d+):([e\d\.\-]+),?//;
		
		$tref->{'distances'}[$id] = $dist;
		push(@{$tref->{'children'}[$vid]}, $id);
	    }
	} else {
	    while ($twig =~ /^(\d+)/) {
		my $id = $1;
		my $dist = 1;

		$twig =~ s/^(\d+),?//;
		
		$tref->{'distances'}[$id] = $dist;
		push(@{$tref->{'children'}[$vid]}, $id);
	    }
	}
	
	&calc_distances($tref, $vid);
    }

    for (my $i = 2; $i <= $#TRANSLATIONS; $i++) {
	for (my $j = 1; $j <= $i - 1; $j++) {
	    $AVERAGES[$i][$j] += $tref->{'distmatrix'}[$i][$j];
	}
    }

    foreach my $a (@{$tref->{'addons'}}) {
	if ($a =~ /^&W\s(.*)$/) {
	    $tref->{'weight'} = $a;
	} elsif ($a =~ /Trees found in bootstrap replicate #\d+/) {
	    $tref->{'bs_source'} = $a;
	}
    }

    &print_tree($tref);
}

sub print_tree {
    my ($tref) = @_;
    
    print $OUT "[" . $tref->{'bs_source'} . "]\n" if $tref->{'bs_source'};
    print $OUT "[! " . $tref->{'name'} . "]\n";
    
    if ($tref->{'weight'}) {
	printf $OUT "begin distances [%s];\n", $tref->{'weight'};
    } else {
        print $OUT "begin distances;\n";
    }
    
    print $OUT "  format triangle=lower nodiagonal nolabels;\n";
    print $OUT "  matrix\n";
    &print_LTmatrix($tref->{'distmatrix'});
    print $OUT "  ;\n";
    print $OUT "end;\n";
}

sub calc_distances {
    my ($tref, $vid) = @_;
    
    my $nochildren = $#{$tref->{'children'}[$vid]};

    for (my $i = 0; $i <= $nochildren; $i++) {
	for (my $j = $i + 1; $j <= $nochildren; $j++) {
	    &calc_distance($tref,
		$tref->{'children'}[$vid][$i], 0.0,
		$tref->{'children'}[$vid][$j], 0.0);
	}
    }
}

sub calc_distance {
    my ($tref, $id1, $dist1, $id2, $dist2) = @_;
    
    if ($id1 > $#TRANSLATIONS && $id2 > $#TRANSLATIONS) {
	foreach my $cid1 (@{$tref->{'children'}[$id1]}) {
	    foreach my $cid2 (@{$tref->{'children'}[$id2]}) {
		&calc_distance($tref,
		    $cid1, $dist1 + $tref->{'distances'}[$id1],
		    $cid2, $dist2 + $tref->{'distances'}[$id2]);
	    }
	}
    } elsif ($id1 > $#TRANSLATIONS) {
	foreach my $cid (@{$tref->{'children'}[$id1]}) {
	    &calc_distance($tref,
		$cid, $dist1 + $tref->{'distances'}[$id1],
		$id2, $dist2);
	}
    } elsif ($id2 > $#TRANSLATIONS) {
	foreach my $cid (@{$tref->{'children'}[$id2]}) {
	    &calc_distance($tref,
		$id1, $dist1,
		$cid, $dist2 + $tref->{'distances'}[$id2]);
	}
    } else {
	if ($id1 > $id2) {
	    die "Internal error: matrix override, aborting ..,\n" if $tref->{'distmatrix'}[$id1][$id2];

	    $tref->{'distmatrix'}[$id1][$id2] =
		$tref->{'distances'}[$id1] + $tref->{'distances'}[$id2] + $dist1 + $dist2;
	} else {
	    die "Internal error: matrix override, aborting ...\n" if $tref->{'distmatrix'}[$id2][$id1];
	    
	    $tref->{'distmatrix'}[$id2][$id1] =
		$tref->{'distances'}[$id1] + $tref->{'distances'}[$id2] + $dist1 + $dist2;
	}
    }
}

##############################################################################
#
# DISTANCE MATRICES --> AVERAGES / STANDARD DEVIATIONS
#
#############################################################################

sub process_dist_file {
    my ($if, $of) = @_;
    
    my $notrees = 0;
    
    print "Processing distances ...\n" if $DEBUG;
    
    $IN = &$OPENFILE($if);
    open($OUT, ">$of");
    
    if ($FILETYPE eq "nexus") {
	print "  Reading taxlabels\n" if $DEBUG;
	
	while (($SECTION ne "distances") && (my $line = &$GETLINE($IN))) {
	    if ($#TRANSLATIONS == -1 && $SECTION eq "taxa") {
		if ($line->{'data'} =~ /^taxlabels/i) {
		    my $i = 1;
		    
		    while ($line->{'data'} =~ /taxlabels\s+('[^']+'|[^\s]+)/i) {
			$TRANSLATIONS[$i++] = $1;
			$line->{'data'} =~ s/taxlabels\s+('[^']+'|[^\s])+(,\s*)?/taxlabels /i;
		    }
		}
	    }
	}
	
	die "No taxlabels found in $if, aborting ...\n" if $#TRANSLATIONS == -1;
    }
    
    for (my $i = 0; $i < (($#TRANSLATIONS ** 2) - $#TRANSLATIONS) / 2; $i++) {
	$AVERAGES[$i] = 0;
	$STDDEVS[$i] = 0;
    }
    
    print "  Calculating averages " if $DEBUG;
    while (my $line = &$GETLINE($IN)) {
	if (($SECTION eq "distances") && ($line->{'data'} =~ /matrix\s+(([e\d\.\-]+\s*)+)/)) {
	    print "." if $DEBUG;
	    my @values = split(/\s+/, $1);
	    my $weight = 1.0;
	    my $bs = 0;
	    
	    if ($OPT_USE_WEIGHTS) {
	        foreach my $a (@{$line->{'addons'}}) {
		   if ($a =~ /^&W\s+([\d\.\-]+)\/([\d\.\-]+)/) {
			$weight = $1 / $2;
		    } elsif ($a =~ /^&W\s+([e\d\.\-]+)/) {
			$weight = $1;
		    } elsif ($a =~ /^Trees found in bootstrap replicate #(\d+)/) {
			$bs = 1;
		    }
		}
	    }
	    
	    $notrees++ if (! $OPT_USE_WEIGHTS) || $weight == 1.0 || $bs;
	    
	    for (my $i = 0; $i <= $#values; $i++) {
		$AVERAGES[$i] += $values[$i] * $weight;
	    }
	}
    }
    print "\n" if $DEBUG;

    die "No distances matrices found in $if, aborting ...\n" if $notrees == -1;
    
    for (my $i = 0; $i <= $#AVERAGES; $i++) {
	$AVERAGES[$i] /= $notrees;
    }
    
    close($IN); $IN = &$OPENFILE($if);

    print "  Calculating standard deviations\n" if $DEBUG;
    while (my $line = &$GETLINE($IN)) {
	if (($SECTION eq "distances") && ($line->{'data'} =~ /matrix\s+(([e\d\.\-]+\s*)+)/)) {
	    print "." if $DEBUG;
	    my @values = split(/\s+/, $1);
	    my $weight = 1.0;
	    
	    if ($OPT_USE_WEIGHTS) {
	        foreach my $a (@{$line->{'addons'}}) {
		    if ($a =~ /^&W\s+([\d\.\-]+)\/([\d\.\-]+)/) {
		        $weight = $1 / $2; last;
		    } elsif ($a =~ /^&W\s+([e\d\.\-]+)/) {
		        $weight = $1; last;
		    }
		}
	    }

	    for (my $i = 0; $i <= $#values; $i++) {
	        $VARIANCES[$i] += (($values[$i] - $AVERAGES[$i]) ** 2.0) * $weight;
	    }
	}
    }
    print "\n" if $DEBUG;

    for (my $i = 0; $i <= $#STDDEVS; $i++) {
	$VARIANCES[$i] = $VARIANCES[$i] / $notrees;
	$VARIANCESINV[$i] = 1.0 / ($VARIANCES[$i] ? $VARIANCES[$i] : 0.000001);
	$STDDEVS[$i] = sqrt($VARIANCES[$i]);
    }

    close($IN);

    &print_averages($if, $notrees);
    
    close($OUT);

    print "done.\n";
}

sub print_averages {
    my ($sf, $notrees) = @_;
    
    print $OUT "#NEXUS\n";
    print $OUT "[! Mean pairwise distances from file " . $sf . "]\n";
    
    &print_taxa();
    
    printf $OUT "begin distances [&N %d];\n", $notrees;
    print  $OUT "  format triangle=lower nodiagonal nolabels;\n";
    print  $OUT "  matrix\n";
    &print_ALTmatrix(\@AVERAGES);
    print  $OUT "  ;\n";
    print  $OUT "end;\n";
    
    printf $OUT "begin inv_variances [&N %d];\n", $notrees;
    print  $OUT "  format triangle=lower nodiagonal nolabels;\n";
    print  $OUT "  matrix\n";
    &print_ALTmatrix(\@VARIANCESINV);
    print  $OUT "  ;\n";
    print  $OUT "end;\n";
 
    printf $OUT "begin deviations [&N %d];\n", $notrees;
    print  $OUT "  format triangle=lower nodiagonal nolabels;\n";
    print  $OUT "  matrix\n";
    &print_ALTmatrix(\@STDDEVS);
    print  $OUT "  ;\n";
    print  $OUT "end;\n";
}

##############################################################################
#
# OUTPUT
#
#############################################################################

sub print_taxa {
    print  $OUT "begin taxa;\n";
    printf $OUT "  dimensions ntax=%d;\n", $#TRANSLATIONS;
    print  $OUT "  taxlabels\n";

    for (my $i = 1; $i <= $#TRANSLATIONS; $i++) {
	print $OUT "    $TRANSLATIONS[$i]\n";
    }

    print  $OUT "    ;\n";
    print  $OUT "end;\n";
}

sub print_LTmatrix {
    my ($mref) = @_;
    
    for (my $i = 2; $i <= $#TRANSLATIONS; $i++) {
	print $OUT "\t";

	for (my $j = 1; $j <= $i - 1; $j++) {
	    printf $OUT "%f\t", $$mref[$i][$j];
	}
	
	print $OUT "\n";
    }
}

sub print_ALTmatrix {
    my ($mref) = @_;
    
    my $row = 0;
    my $i = 0;

    foreach my $val (@{$mref}) {
	printf $OUT "\t%f", $val;
	
	if ($i == $row) {
	    print $OUT "\n";

	    $row++;
	    $i = 0;
	} else {
	    $i++;
	}
    }
}

##############################################################################
#
# MAIN
#
#############################################################################

sub main {
    my $if = '';
    my $of_dist = '';
    my $of_avg = '';
    my $do_dist = 0;
    my $do_avg = 0;
    my $tree_first = 1;
    my $tree_last = 0;

    while ($#ARGV > -1) {
	my $arg = shift(@ARGV);

	if ($arg eq "--dist-output") {
	    $of_dist = shift(@ARGV);
	} elsif ($arg eq "--avg-output") {
	    $of_avg = shift(@ARGV);
	} elsif ($arg eq "-f") {
	    $tree_first = shift(@ARGV);
	} elsif ($arg eq "-l") {
	    $tree_last = shift(@ARGV);
	} elsif ($arg eq "-d") {
	    $do_dist = 1;
	} elsif ($arg eq "-a") {
	    $do_avg = 1;
	} elsif ($arg eq "-1") {
	    $OPT_ONE = 1;
	} elsif ($arg eq "-w") {
	    $OPT_USE_WEIGHTS = 1;
	} else {
	    $if = $arg;
	}
    }
    
    die "Usage: $0 [option ...] (infile|-)\n" .
	"  Options:\n" .
	"  --dist-output file   Use file for distance matrices output\n" .
	"  --avg-output file    Use file for averages output\n" .
	"  -f no                Only use trees starting from no\n" .
	"                       (first tree is numbered 1)\n" .
	"  -t no                Only use trees up to no\n" .
	"  -d                   Calculate distance matrices\n" .
	"  -a                   Calculate averages / standard deviations\n" .
	"  -1                   Set all distances to 1.0\n" .
	"  -w                   Use bootstrapping weights if found\n" .
	"\n" .
	"  Defaults:\n" .
	"  Both distance matrices and averages / std.dev's are calculated\n" .
	"  All trees found are used\n" .
	"  Bootstrapping weights are ignored\n"
	unless $if;
    
    
    &det_filetype($if);
    
    print "$if seems to be a $FILETYPE file.\n";

    if (! ($do_dist || $do_avg)) {
	$do_dist = 1; $do_avg = 1;
    }

    my $base = 'pardist';
    if ($if =~ /^(.*?)(_(distances|averages))?\.[^\.]+$/) {
	$base = $1;
    }
    
    if ($do_dist) {
	if ($FILETYPE eq 'nexus') {
	    $of_dist = $base . "_distances.nex" unless $of_dist;
	    
	    &process_tree_file($if, $of_dist, $tree_first, $tree_last);
	    
	    $if = $of_dist; &det_filetype($if);
	} elsif ($FILETYPE eq 'phylip') {
	    &iop_phylip2nexus($if);
	}
    }
    
    if ($do_avg) {
	$of_avg = $base . "_averages.nex" unless $of_avg;
	
	&process_dist_file($if, $of_avg);
    }
    
    return 0;
}

exit(&main());


