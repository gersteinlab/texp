#!/usr/bin/perl
#####################################
# Program: remove_homopol.pl  -  Date: Fri Nov  7 13:45:24 EST 2014
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

while ( <> ) {
	chomp $_;
	if ( $_ =~ /^[@]/ ) { 
		print "$_\n";
	} else {	
		my @tokens = split ( /[ \t\cI]+/,$_ );
	
		my %nts_count = ();
		my @nts = split(//,$tokens[9]);
		foreach my $nt (@nts) {
			$nts_count{$nt}++;
		}
		my $sum = 0;
		foreach my $nt ( keys %nts_count ) {
			$sum += $nts_count{$nt};
		}
		my $print = 1;
		foreach my $nt ( keys %nts_count ) {
			if ( $nts_count{$nt}/$sum >= 0.70 ) {
				$print = 0;
			}
		}
		if ( $print ) {
			print "$_\n";
		}
	}
}
