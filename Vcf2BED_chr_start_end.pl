#!/usr/bin/perl -i
#Vcf2BED_chr_start_end.pl
#by HIRAO Akira

while ($line = <>) {
	chomp $line;
	if($line !~ m/^#/){	    
		($Chr, $end, @info) = split /\s+/, $line; 
	    $start = $end  - 1;
		print $Chr, "\t", $start, "\t", $end, "\n";
	}
}
