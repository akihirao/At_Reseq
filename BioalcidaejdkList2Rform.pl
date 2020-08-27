#!/usr/bin/perl -i
#BioalcidaejdkList2Rform.pl
#by HIRAO Akira


#当該SNVのポジションを
while ($line = <>) {
	chomp $line;
	($Chr, $position, $ref, $sample, $ref_dash, $alt) = split /\s+/, $line; 
#	$ref = substr($ref_orig, -1, 1, );
#	$alt = substr($alt_dash, -1, 1, );

	chop $ref;
	chop $alt;

	$mutation_length = length($alt) - length($ref);
	if($mutation_length > 0){
		$mutation_type = "Insertion";
	}elsif($mutation_length < 0){
		$mutation_type = "Deletion";
	}else{
		$mutation_type = "SNP"; 
	}

	$start_position = $end_position -1;
	print $Chr, "\t", $position, "\t", $ref, "\t", $alt, "\t", $mutation_type, "\t", $mutation_length, "\t", $sample, "\n";
	
}#close while()
