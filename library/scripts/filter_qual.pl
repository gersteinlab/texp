#!/usr/bin/perl
#####################################
# Program: filter_qual.pl  -  Date: Thu Nov 20 13:23:15 EST 2014
# Autor: Fabio C. P. Navarro - Ludwig
# Goal:
#
# Input:
#
# Output:
#
#####################################


use Getopt::Long;
use strict;
use warnings;

#Remove alignments with low alignment score and remove alignments with clipping on both sides.
#GA-DELTA1619-1_0006:8:9:17662:17830#0/1 0       chr1    11057   0       7S22M46S        *       0       0       CCCTGTCCCGGGCTGGGGCGGGGGGAGGGGAAGCTTGTTGGGTGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG       \HFRU_SOOSPRLUBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB     AS:i:44 XS:i:44 XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:22 YT:Z:UU
#open ( IN0, "<$file" );
while ( <> ) {
	if ( $_ =~ /^[@]/ ) {
		print $_;
	} else { 
		chomp $_;
		my @tokens = split ( /[ \t\cI]+/,$_ );
	
		#Deal with CIGAR
		my $cigar = $tokens[5];
		my $num = "";
		my $count = 0;
		my @tmp = ();
	
		my $remove_clipping = 0;
	
		foreach my $char ( split(//,$tokens[5]) ) {
			if ( $char !~ /[0-9]/ ) {
				$tmp[$count]->[0] = int($num);
				$tmp[$count]->[1] = $char;
				$count++;
				$num = "";
			} else {
				$num .= $char;
			}
	 	}
		
		if ( $tmp[0]->[1] eq "S" && $tmp[0]->[0] > 2 && $tmp[-1]->[1] eq "S" && $tmp[-1]->[0] > 2 ) {
			$remove_clipping = 1;
		} else {
			$remove_clipping = 0;
		}
	
		#Deal with alignment scores. 
		my $alignment_score_threshold = length($tokens[9])*2*0.80;
	
		my $remove_alignment_score = 0;
	
		for ( my $i = 11; $i < scalar(@tokens); $i++ ) {
			if ( $tokens[$i] =~ /^AS:i:/ ) {
				my @score = split(/[\:]/,$tokens[$i]);
				if ( $score[2] < $alignment_score_threshold ) {
					$remove_alignment_score = 1;
				}
				else {
					$remove_alignment_score = 0;
				}
				last;
			}
		}
	
		if ( !$remove_alignment_score && !$remove_clipping ) {
			print "$_\n";
		}
	}
}
