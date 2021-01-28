#!/usr/bin/perl -i
#VcfHomoHeteroCount.pl
#by HIRAO Akira

#list of SampleID 
@sample = ("A011","A012","A014","A021","A023","A024","A031","A033","A034","A113","A114","A115","A121","A123","A125","A131","A132","A135","A212","A213","A215","A221","A222","A225","A231","A233","A234","A311","A313","A314","A321","A323","A324","A331","A332","A334");


$NoSample = @sample;

open(OUT, ">SnpEff.AT.gene.unsorted.txt");
print OUT "Chr","\t","Position","\t","Ref","\t","Alt","\t","Gene","\t","Annotation","\t","Annotation_impact","\t","Annotation_gene_name","\t","Annotation_full","\n";

for($i = 0; $i < $NoSample; $i++){

	$target_sample = $sample[$i];

	#on Takeru
	$source_file_path = "/zfs/Arabidopsis/work/At_Reseq/vcf_out/".$target_sample."/".$target_sample.".final.mutations.snpeff.vcf";

	open(SOURCE, $source_file_path);

	$NoHetero = 0; 
	$NoHomo = 0;
	$NoTotalMutation = 0;

	while ($line = <SOURCE>) {
		chomp $line;
		if($line !~ m/^#/){
			$NoTotalMutation = $NoTotalMutation + 1;

			($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SampleInfo) = split /\s+/, $line; 
		
			@INFO_vec = split /;/, $INFO;

			for($count = 0; $count < @INFO_vec; $count++){
				
				$target_INFO = $INFO_vec[$count];
				if($target_INFO =~ m/^ANN/){

					@SnpEff_full_info_vec = split /\|/, $target_INFO;
					$SnpEff_Annotation = $SnpEff_full_info_vec[1];
					$SnpEff_Annotation_Impact = $SnpEff_full_info_vec[2];
					$SnpEff_Annotation_Gene_Name = $SnpEff_full_info_vec[3];
					$SnpEff_Annotation_Gene_ID = $SnpEff_full_info_vec[4];

					#($GT, $AD, $DP, $GQ, $PGT, $PID, $PL, $PS) = split /:/, $SampleInfo;
					print OUT $CHROM, "\t";
					print OUT $POS, "\t";
					print OUT $REF, "\t";
					print OUT $ALT, "\t";
					if($SnpEff_Annotation_Impact eq "MODIFIER"){
						print OUT "\t", "\t";
					}else{
						print OUT $SnpEff_Annotation_Gene_ID, "\t";
						print OUT $SnpEff_Annotation, "\t";
					}
					print OUT $SnpEff_Annotation_Impact, "\t";					
					print OUT $SnpEff_Annotation_Gene_Name, "\t";
					print OUT $SnpEff_Annotation, "\n";
				}

			}			
		}
	}
}

close(OUT);

