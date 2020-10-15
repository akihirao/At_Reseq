#!/usr/bin/perl -i
#Vcf2Position.pl
#by HIRAO Akira

while ($line = <>) {
	chomp $line;
	if($line !~ m/^#/){	    
		($Chr, $end, @info) = split /\s+/, $line; 
	    $start = $end  - 1;
		print $Chr, ":",$end, "\n";
	}
}
