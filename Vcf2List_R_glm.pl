#!/usr/bin/perl -i
#Vcf2List_R_glm.pl
#by HIRAO Akira

#list of SampleID 
@sample = ("A011","A012","A014","A021","A023","A024","A031","A033","A034","A113","A114","A115","A121","A123","A125","A131","A132","A135","A212","A213","A215","A221","A222","A225","A231","A233","A234","A311","A313","A314","A321","A323","A324","A331","A332","A334");

@dose = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0);
@treatment = ("Control","Control","Control","Control","Control","Control","Control","Control","Control","Low","Low","Low","Low","Low","Low","Low","Low","Low","Middle","Middle","Middle","Middle","Middle","Middle","Middle","Middle","Middle","High","High","High","High","High","High","High","High","High");

$NoSample = @sample;


print "Chr", "\t", "Position", "\t", "Ref", "\t", "Alt", "\t", "Type", "\t", "Length", "\t", "Sample1", "\t", "Zygosity1", "\t", "Sample2", "\t", "Zygosity2", "\t", "Sample3", "\t", "Zygosity3", "\t", "Treatment", "\t", "Dose", "\n";


while ($line = <>) {
	chomp $line;
	if($line !~ m/^#/){

		($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @SampleInfo) = split /\s+/, $line; 
		$Type = "";
		
		$REF_len = length($REF); $ALT_len = length($ALT);
		$Length = $ALT_len - $REF_len;

		#categorizing mutation types; SBS, deletion, insertion
		 if($REF_len == 1 && $ALT_len == 1){
		 	$Type = "SBS";
		 }elsif($REF_len > $ALT_len){
		 	$Type = "Deletion";
		 }else{
		 	$Type = "Insertion";
		 }

		print $CHROM, "\t", $POS, "\t", $REF, "\t", $ALT, "\t", $Type, "\t", $Length, "\t";

		($A011_info,$A012_info,$A014_info,$A021_info,$A023_info,$A024_info,$A031_info,$A033_info,$A034_info,$A113_info,$A114_info,$A115_info,$A121_info,$A123_info,$A125_info,$A131_info,$A132_info,$A135_info,$A212_info,$A213_info,$A215_info,$A221_info,$A222_info,$A225_info,$A231_info,$A233_info,$A234_info,$A311_info,$A313_info,$A314_info,$A321_info,$A323_info,$A324_info,$A331_info,$A332_info,$A334_info) = @SampleInfo;

		@A011_GT_info = split /:/, $A011_info;
		$A011_GT_raw = $A011_GT_info[0];
		$A011_GT = substr($A011_GT_raw,0,1) + substr($A011_GT_raw,2,1);

		@A012_GT_info = split /:/, $A012_info;
		$A012_GT_raw = $A012_GT_info[0];		
		$A012_GT = substr ($A012_GT_raw,0,1) + substr ($A012_GT_raw,2,1);

		@A014_GT_info = split /:/, $A014_info;
		$A014_GT_raw = $A014_GT_info[0];
		$A014_GT = substr ($A014_GT_raw,0,1) + substr ($A014_GT_raw,2,1);

		@A021_GT_info = split /:/, $A021_info;
		$A021_GT_raw = $A021_GT_info[0];
		$A021_GT = substr ($A021_GT_raw,0,1) + substr ($A021_GT_raw,2,1);

		@A023_GT_info = split /:/, $A023_info;
		$A023_GT_raw = $A023_GT_info[0];
		$A023_GT = substr ($A023_GT_raw,0,1) + substr ($A023_GT_raw,2,1);

		@A024_GT_info = split /:/, $A024_info;
		$A024_GT_raw = $A024_GT_info[0];
		$A024_GT = substr ($A024_GT_raw,0,1) + substr ($A024_GT_raw,2,1);

		@A031_GT_info = split /:/, $A031_info;
		$A031_GT_raw = $A031_GT_info[0];
		$A031_GT = substr ($A031_GT_raw,0,1) + substr ($A031_GT_raw,2,1);

		@A033_GT_info = split /:/, $A033_info;
		$A033_GT_raw = $A033_GT_info[0];
		$A033_GT = substr ($A033_GT_raw,0,1) + substr ($A033_GT_raw,2,1);

		@A034_GT_info = split /:/, $A034_info;
		$A034_GT_raw = $A034_GT_info[0];
		$A034_GT = substr ($A034_GT_raw,0,1) + substr ($A034_GT_raw,2,1);

		@A113_GT_info = split /:/, $A113_info;
		$A113_GT_raw = $A113_GT_info[0];
		$A113_GT = substr ($A113_GT_raw,0,1) + substr ($A113_GT_raw,2,1);

		@A114_GT_info = split /:/, $A114_info;
		$A114_GT_raw = $A114_GT_info[0];
		$A114_GT = substr ($A114_GT_raw,0,1) + substr ($A114_GT_raw,2,1);

		@A115_GT_info = split /:/, $A115_info;
		$A115_GT_raw = $A115_GT_info[0];
		$A115_GT = substr ($A115_GT_raw,0,1) + substr ($A115_GT_raw,2,1);

		@A121_GT_info = split /:/, $A121_info;
		$A121_GT_raw = $A121_GT_info[0];
		$A121_GT = substr ($A121_GT_raw,0,1) + substr ($A121_GT_raw,2,1);

		@A123_GT_info = split /:/, $A123_info;
		$A123_GT_raw = $A123_GT_info[0];
		$A123_GT = substr ($A123_GT_raw,0,1) + substr ($A123_GT_raw,2,1);

		@A125_GT_info = split /:/, $A125_info;
		$A125_GT_raw = $A125_GT_info[0];
		$A125_GT = substr ($A125_GT_raw,0,1) + substr ($A125_GT_raw,2,1);

		@A131_GT_info = split /:/, $A131_info;
		$A131_GT_raw = $A131_GT_info[0];
		$A131_GT = substr ($A131_GT_raw,0,1) + substr ($A131_GT_raw,2,1);

		@A132_GT_info = split /:/, $A132_info;
		$A132_GT_raw = $A132_GT_info[0];
		$A132_GT = substr ($A132_GT_raw,0,1) + substr ($A132_GT_raw,2,1);

		@A135_GT_info = split /:/, $A135_info;
		$A135_GT_raw = $A135_GT_info[0];
		$A135_GT = substr ($A135_GT_raw,0,1) + substr ($A135_GT_raw,2,1);

		@A212_GT_info = split /:/, $A212_info;
		$A212_GT_raw = $A212_GT_info[0];
		$A212_GT = substr ($A212_GT_raw,0,1) + substr ($A212_GT_raw,2,1);

		@A213_GT_info = split /:/, $A213_info;
		$A213_GT_raw = $A213_GT_info[0];
		$A213_GT = substr ($A213_GT_raw,0,1) + substr ($A213_GT_raw,2,1);

		@A215_GT_info = split /:/, $A215_info;
		$A215_GT_raw = $A215_GT_info[0];
		$A215_GT = substr ($A215_GT_raw,0,1) + substr ($A215_GT_raw,2,1);

		@A221_GT_info = split /:/, $A221_info;
		$A221_GT_raw = $A221_GT_info[0];
		$A221_GT = substr ($A221_GT_raw,0,1) + substr ($A221_GT_raw,2,1);

		@A222_GT_info = split /:/, $A222_info;
		$A222_GT_raw = $A222_GT_info[0];
		$A222_GT = substr ($A222_GT_raw,0,1) + substr ($A222_GT_raw,2,1);

		@A225_GT_info= split /:/, $A225_info;
		$A225_GT_raw = $A225_GT_info[0];
		$A225_GT = substr ($A225_GT_raw,0,1) + substr ($A225_GT_raw,2,1);

		@A231_GT_info = split /:/, $A231_info;
		$A231_GT_raw = $A231_GT_info[0];
		$A231_GT = substr ($A231_GT_raw,0,1) + substr ($A231_GT_raw,2,1);

		@A233_GT_info = split /:/, $A233_info;
		$A233_GT_raw = $A233_GT_info[0];
		$A233_GT = substr ($A233_GT_raw,0,1) + substr ($A233_GT_raw,2,1);

		@A234_GT_info = split /:/, $A234_info;
		$A234_GT_raw = $A234_GT_info[0];
		$A234_GT = substr ($A234_GT_raw,0,1) + substr ($A234_GT_raw,2,1);

		@A311_GT_info = split /:/, $A311_info;
		$A311_GT_raw = $A311_GT_info[0];
		$A311_GT = substr ($A311_GT_raw,0,1) + substr ($A311_GT_raw,2,1);

		@A313_GT_info = split /:/, $A313_info;
		$A313_GT_raw = $A313_GT_info[0];
		$A313_GT = substr ($A313_GT_raw,0,1) + substr ($A313_GT_raw,2,1);

		@A314_GT_info = split /:/, $A314_info;
		$A314_GT_raw = $A314_GT_info[0];
		$A314_GT = substr ($A314_GT_raw,0,1) + substr ($A314_GT_raw,2,1);

		@A321_GT_info = split /:/, $A321_info;
		$A321_GT_raw = $A321_GT_info[0];
		$A321_GT = substr ($A321_GT_raw,0,1) + substr ($A321_GT_raw,2,1);

		@A323_GT_info = split /:/, $A323_info;
		$A323_GT_raw = $A323_GT_info[0];
		$A323_GT = substr ($A323_GT_raw,0,1) + substr ($A323_GT_raw,2,1);

		@A324_GT_info = split /:/, $A324_info;
		$A324_GT_raw = $A324_GT_info[0];
		$A324_GT = substr ($A324_GT_raw,0,1) + substr ($A324_GT_raw,2,1);

		@A331_GT_info = split /:/, $A331_info;
		$A331_GT_raw = $A331_GT_info[0];
		$A331_GT = substr ($A331_GT_raw,0,1) + substr ($A331_GT_raw,2,1);

		@A332_GT_info = split /:/, $A332_info;
		$A332_GT_raw = $A332_GT_info[0];
		$A332_GT = substr ($A332_GT_raw,0,1) + substr ($A332_GT_raw,2,1);

		@A334_GT_info = split /:/, $A334_info;
		$A334_GT_raw = $A334_GT_info[0];
		$A334_GT = substr ($A334_GT_raw,0,1) + substr ($A334_GT_raw,2,1);
#			print $A334_GT, "\t";


		@GT_list = ($A011_GT,$A012_GT,$A014_GT,$A021_GT,$A023_GT,$A024_GT,$A031_GT,$A033_GT,$A034_GT,$A113_GT,$A114_GT,$A115_GT,$A121_GT,$A123_GT,$A125_GT,$A131_GT,$A132_GT,$A135_GT,$A212_GT,$A213_GT,$A215_GT,$A221_GT,$A222_GT,$A225_GT,$A231_GT,$A233_GT,$A234_GT,$A311_GT,$A313_GT,$A314_GT,$A321_GT,$A323_GT,$A324_GT,$A331_GT,$A332_GT,$A334_GT);

		$last_sample_ID = ();
		$sample_out_count = 0;
		for($i = 0; $i < $NoSample; $i++){
			$target_GT = $GT_list[$i];
			if($target_GT > 0){
				$target_out_sample = $sample[$i];
				print $target_out_sample, "\t";
				if($target_GT == 1){
					print "hetero", "\t";
				}else{
					print "homo", "\t";
				}
			$sample_out_count++;
			$last_sample_ID = $i;
			}
		}
		#}
		#foreach $out_gt (@GT_list){
		#	print $out_gt, "\t";
		#}

		if($sample_out_count == 1){
			print "NA", "\t", "NA", "\t", "NA", "\t", "NA", "\t";
		}elsif(($sample_out_count == 2)){
			print "NA", "\t","NA", "\t";
		}else{

		}

		$dose_output = $dose[$last_sample_ID];
		$treatment_output = $treatment[$last_sample_ID];

		print $treatment_output, "\t", $dose_output, "\n";

	}
}






