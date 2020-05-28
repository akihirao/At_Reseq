#!/usr/bin/perl -i
#FilteredText2BED.pl
#from 2019,4.8 to 2019,5.22
#by HIRAO Akira


#当該SNVのポジションを
while ($line = <>) {
	chomp $line;
	($Chr, $end_position, @info) = split /\s+/, $line; 
	$start_position = $end_position -1;
	print $Chr, "\t", $start_position, "\t", "$end_position", "\n";
	
}#close while()
