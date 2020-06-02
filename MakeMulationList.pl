#!/usr/bin/perl -i
#MakeMulationList.pl
#by HIRAO Akira

if(@ARGV <2){
	die "USAGE: perl MakeMutationList.pl input.vcf sample_name\n";
}

($input_vcf, $sample_name) = @ARGV;

open (SOURCE, $input_vcf);

while ($line = <SOURCE>) {
	chomp $line;
	if($line =~ m/^#/){
		#print $line, "\n";
	}else{		
		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SampleInfo) = split /\s+/, $line; 
		
		$no_REF = length($REF);
		$no_ALT = length($ALT);
		$mutation_bp = $no_ALT - $no_REF;
	
		if ($no_REF == 1 && $no_ALT > 1){
			$mutation_type = "Insertion";
		}elsif($no_REF > 1 && $no_ALT == 1){
			$mutation_type = "Deletion";
		}else{
			$mutation_type = "SNP";
		}
				
		($GT, $AD, $DP, $GQ, $PGT, $PID, $PL, $PS) = split /:/, $SampleInfo;
	
		print $CHROM, "\t", $POS, "\t", $REF, "\t", $ALT, "\t", $mutation_type, "\t", $mutation_bp, "\t", $sample_name, "\t", $GQ, "\n";		
			
	}
}

close (SOURCE);