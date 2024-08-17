#!/usr/bin/perl
use strict;
use warnings;

use Mods::GenoMetaAss qw(systemW gzipwrite gzipopen);

#####################
#This script takes contig coverage file, rounds the numbers to the nearest intiger, and use that information to duplicate the contigs in fasta file; In the second part it takes extended contigs and split them based on desiger fragment size:
# Check if provided input and output filenames as command-line arguments
#Run with: ./split_fasta4metaMDBG.pl scaffolds.fasta.filt Coverage.percontig scaffolds.fasta.extendcont_temp 150 scaffolds.fasta.extendcont.split
#####################

#v0.1: 17.3.24 (FH): changed to fq.gz output

if (@ARGV <= 2) {
	die "Usage: $0 <fasta_file> <key_value_file> <output_fasta_file_temp> <fragment_length> <output_fastq_file>\n";
}

my $fasta_file = $ARGV[0]; #scaffolds.fasta.filt
my $key_value_file = $ARGV[1]; #Coverage.percontig
#my $output_fasta_file_temp = $ARGV[2]; #scaffolds.fasta.extendcont_temp
#my $fragment_length = $ARGV[3]; #fragment length, fragment size based on which you want to split your extended contigs
my $output_fasta_file = $ARGV[2]; #scaffolds.fasta.extendcont.split

#define constants
my $maxCtgL = 20000;
my $SplitCtgL = 15000;




#######################
#First part of script:  Multiplication of contigs
#######################

print "Subsplitting assembly in prep for metaMDBG\nInput: $fasta_file, $key_value_file\nOutput: $output_fasta_file\n";

# Read key-value pairs from the input file
my %key_values;
#open my $key_value_fh, '<', $key_value_file or die "Cannot open file '$key_value_file' for reading: $!\n";
my ($key_value_fh,$OK) = gzipopen($key_value_file,"contig abundances",1);


# Read the input file line by line
while (<$key_value_fh>) {
	chomp; # Remove newline characters

# Split the line into columns
	my ($key, $value) = split /\s+/, $_, 2;

# Round the value to the nearest integer, and modify key by adding >
	#my $rounded_value = int($value + 0.5);
	#my $modified_key = ">$key";

# Assign the key and rounded value to the hash
	$key_values{"$key"} = $value;
}

# Check the pairs of keys and values:
#foreach my $key (keys %key_values) {
#	print "$key\t$key_values{$key}\n";
#}
#or
#print "@{[%key_values]}";


# Open the input FASTA file for reading
open my $fasta_fh, '<', $fasta_file or die "Cannot open file '$fasta_file' for reading: $!\n";

# Open the output FASTA file for writing
#open my $output1_fh, '>', $output_fasta_file or die "Cannot open file '$output_fasta_file' for writing: $!\n";
my $output1_fh = gzipwrite($output_fasta_file,"split4metaMDBG fq out");


# Initialize variables to store header and sequence
my $header = "";
my $sequence = "";

my $covL = 0;

# Read the input file line by line
while (<$fasta_fh>) {
	chomp; # Remove newline characters

# Check if the line is a header line
	if (/^>/) {
		# If the header is not empty, print the previous sequence and check that keys are matching
		if ($header ne "" ){
			if (exists $key_values{$header}) { 
				$covL = $key_values{$header};
			} else {;}#print STDERR "could not find coverage for $header\n";}
			SPLT_sequence($header, $sequence,$covL,$output1_fh);
		}
		# Update the header and reset the sequence
		$header = substr($_,1);
		$sequence = "";
		$covL=0;
		next;
	} 
	# Concatenate sequence lines
	$sequence .= $_;
}

#Process the last sequence in the file
if ($header ne "" ) { 
	SPLT_sequence($header, $sequence,$covL,$output1_fh);
}

# Close the filehandles
close $fasta_fh;
close $key_value_fh;
close $output1_fh;


print "Done resplitting fasta for metaMDBG to file $output_fasta_file\n\n";

exit(0);









# Function to process a sequence, split it into fragments, and write to the output file
sub SPLT_sequence {
	my ($header2, $sequence2, $coverage, $output2_fh) = @_; #fragment length is 150bp for this trial
	
	return if ($coverage <= 0);
	
	my $tot_length = length($sequence2);

# Calculate the number of fragments
	my $num_fragments = 1;
	if ($tot_length > $maxCtgL){
		$num_fragments = int($tot_length/$SplitCtgL + 0.5);
	}
	my $fragment_length = $tot_length/$num_fragments;
	
	
	my $numCps = int($coverage+0.5);
	
	
	#print "$header2\t$tot_length\t$num_fragments\t$fragment_length\t$numCps\n";

# Split the sequence into fragments and write to the output file
	for (my $cp=0;$cp < $numCps; $cp++){
		for my $i (0 .. ($num_fragments -1) ) {
			my $start = $i * $fragment_length;
			my $end = ($i + 1) * $fragment_length - 1;
			if ($i == ($num_fragments -1)){
				$end = $tot_length;
			}
			
			if ($end-$start < 800 && ($num_fragments -1) != $i){
				print "WAY too short:: $end - $start $header2\n";
			}

			my $fragment_header = "${header2}_CPY${cp}_SPL_$i"; #update the header
			my $fragment_sequence = substr($sequence2, $start, $fragment_length); #new sequence
#FASTA out
#			print $output2_fh "\@$fragment_header\n$fragment_sequence\n";
#FASTQ out
			my $QUALstr = "O" x length($fragment_sequence); #qual ~ 40
			print $output2_fh "\@$fragment_header\n$fragment_sequence\n+\n$QUALstr\n";
		}

	}
	#sometimes, final fragments can be really short:
	# If the final fragment is shorter than $fragment_length, add it to the previous one (doesn't work; yes it still puts it in the line, but it doesn't add it to previous one....)
	#if ($num_fragments * $fragment_length < $tot_length) {
	#	my $final_fragment_header = "$header2\_SPL_$num_fragments"; 
	#	my $final_fragment_sequence = substr($sequence2, $num_fragments * $fragment_length);

	#	print $output2_fh "$final_fragment_header\n$final_fragment_sequence\n";
	#}
}






