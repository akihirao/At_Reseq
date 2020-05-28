#!/usr/bin/perl -i
#FilteredNoCallVcf.pl
#from 2019,4.8 to 2019,5.22
#by HIRAO Akira


#当該SNVのポジションを
while ($line = <>) {
	chomp $line;
	if($line =~ m/^#/){
		print $line, "\n";
	}else{		
		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SampleInfo) = split /\s+/, $line; 
		if($SampleInfo =~ m/^0/){
			print $line, "\n";
		}else{

		}
	}
}#close while()
