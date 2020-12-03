#!/usr/bin/perl -i
#VcfHomoHeteroCount.pl
#by HIRAO Akira

#list of SampleID 
@sample = ("A011","A012","A014","A021","A023","A024","A031","A033","A034","A113","A114","A115","A121","A123","A125","A131","A132","A135","A212","A213","A215","A221","A222","A225","A231","A233","A234","A311","A313","A314","A321","A323","A324","A331","A332","A334");


$NoSample = @sample;

open(OUT, ">AT.M2.homo.hetero.count.txt");
print OUT "Sample","\t","No.Mutation", "\t","No.Homo.Mutation", "\t","No.Hetero.Mutation","\n";

for($i = 0; $i < $NoSample; $i++){

	$target_sample = $sample[$i];

	#on Takeru
	$source_file_path = "/zfs/Arabidopsis/work/At_Reseq/vcf_out/".$target_sample."/".$target_sample.".homo.hetero.familyclustered.vcf";
#	$source_file_path = $target_sample."/".$target_sample.".homo.hetero.familyclustered.vcf";

	open(SOURCE, $source_file_path);

	$NoHetero = 0; 
	$NoHomo = 0;
	$NoTotalMutation = 0;

	while ($line = <SOURCE>) {
		chomp $line;
		if($line !~ m/^#/){
			$NoTotalMutation = $NoTotalMutation + 1;

			($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SampleInfo) = split /\s+/, $line; 
		
			($GT, $AD, $DP, $GQ, $PGT, $PID, $PL, $PS) = split /:/, $SampleInfo;
		
			if(($GT eq "1/1")||($GT eq "1|1")){
				$NoHomo = $NoHomo + 1;
			}else{
				$NoHetero = $NoHetero + 1;
			}
		}
	}

	print OUT $target_sample, "\t", $NoTotalMutation, "\t", $NoHomo, "\t", "$NoHetero", "\n";
#	print $NoTotalMutation, "\t", $NoHomo, "\t", "$NoHetero", "\n";

}

close(OUT);

