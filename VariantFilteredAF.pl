#!/usr/bin/perl -i
#VariantFilteredAF.pl
#by HIRAO Akira


#definign threshold value for proportion of mutant reads at a site
$AF_threshold = 0.25;

while ($line = <>) {
	chomp $line;
	if($line =~ m/^#/){
		print $line, "\n";
	}else{		
		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SampleInfo) = split /\s+/, $line; 
		
		($GT, $AD, $DP, $GQ, $PGT, $PID, $PL, $PS) = split /:/, $SampleInfo;
		($AD0, $AD1) = split /,/, $AD;
		$AF = $AD1 / ($AD0 + $AD1); 	
		
		if($AF > $AF_threshold){
			print $line, "\n";		
		}else{

		}
	}
}
