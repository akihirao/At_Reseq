#!/usr/bin/perl -i
#BioalcidaejdkOut2BED.pl
#by HIRAO Akira


#当該SNVのポジションを
while ($line = <>) {
	chomp $line;
	($Chr, $end_position, @info) = split /\s+/, $line; 
	$start_position = $end_position -1;
	print $Chr, "\t", $start_position, "\t", "$end_position", "\n";
	
}#close while()
