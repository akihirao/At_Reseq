#!/usr/bin/perl -i
#VariantFilteredAF_homo.pl
#selecting homozygous mutations whrere proportion of mutant reads > 80%
#by HIRAO Akira


#definign threshold value for proportion of mutant reads at a site
$AF_low_threshold = 0.25;
$AF_high_threshold = 0.8;

while ($line = <>) {
	chomp $line;
	if($line =~ m/^#/){
		print $line, "\n";
	}else{		
		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SampleInfo) = split /\s+/, $line; 
		
		($GT, $AD, $DP, $GQ, $PGT, $PID, $PL, $PS) = split /:/, $SampleInfo;
		if($GT eq "./."){
		}else{	
			($AD0, $AD1) = split /,/, $AD;
			$AF = $AD1 / ($AD0 + $AD1); 	
			
			if($AF >= $AF_high_threshold){
				print $line, "\n";		
			}

		}
	}
}
