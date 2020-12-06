#!/usr/bin/perl -i
#CoverageSummarizingMosdepth.pl
#by HIRAO Akira


#Total bases of TAIR10
$total_base = 119668634; 

$Chr1_len = 30427671;
$Chr2_len = 19698289;
$Chr3_len = 23459830;
$Chr4_len = 18585056;
$Chr5_len = 26975502;

$total_nuclear_chr_len = $Chr1_len + $Chr2_len + $Chr3_len + $Chr4_len + $Chr5_len;


#list of SampleID 
@sample = ("AT01","AT02","AT03","AT11","AT12","AT13","AT21","AT22","AT23","AT31","AT32","AT33","A011","A012","A014","A021","A023","A024","A031","A033","A034","A113","A114","A115","A121","A123","A125","A131","A132","A135","A212","A213","A215","A221","A222","A225","A231","A233","A234","A311","A313","A314","A321","A323","A324","A331","A332","A334");



$NoSample = @sample;

open(OUT, ">AT48.coverage.summary.txt");
print OUT "Sample","\t","Coverage.total", "\t","Coverage.nDNA","\t","Genome.Coverage.rate.x10.total","\n";


#for($i = 0; $i < 2; $i++){
for($i = 0; $i < $NoSample; $i++){

	$target_sample = $sample[$i];

	$source_file1_path = "/zfs/Arabidopsis/work/At_Reseq/bwa_out/".$target_sample."/".$target_sample.".mosdepth.summary.txt";
	open(SOURCE1, $source_file1_path);

	$total_nDNA_bases = 0;
	$total_nDNA_length = 0;

	$line = <SOURCE1>; #header 

	#Chr1
	$line = <SOURCE1>; 
	chomp $line;
	($chrom, $length, $bases, $mean, $min, $max) = split /\s+/, $line; 
	$total_nDNA_length = $total_nDNA_length + $length;
	$total_nDNA_bases = $total_nDNA_bases + $bases;

	#Chr2
	$line = <SOURCE1>; 
	chomp $line;
	($chrom, $length, $bases, $mean, $min, $max) = split /\s+/, $line; 
	$total_nDNA_length = $total_nDNA_length + $length;
	$total_nDNA_bases = $total_nDNA_bases + $bases;

	#Chr3
	$line = <SOURCE1>; 
	chomp $line;
	($chrom, $length, $bases, $mean, $min, $max) = split /\s+/, $line; 
	$total_nDNA_length = $total_nDNA_length + $length;
	$total_nDNA_bases = $total_nDNA_bases + $bases;

	#Chr4
	$line = <SOURCE1>; 
	chomp $line;
	($chrom, $length, $bases, $mean, $min, $max) = split /\s+/, $line; 
	$total_nDNA_length = $total_nDNA_length + $length;
	$total_nDNA_bases = $total_nDNA_bases + $bases;
	
	#Chr5
	$line = <SOURCE1>; 
	chomp $line;
	($chrom, $length, $bases, $mean, $min, $max) = split /\s+/, $line; 
	$total_nDNA_length = $total_nDNA_length + $length;
	$total_nDNA_bases = $total_nDNA_bases + $bases;

	$total_nDNA_coverage = $total_nDNA_bases/$total_nDNA_length;
	$total_nDNA_coverage_out = sprintf("%.2f", $total_nDNA_coverage);

	#ChrM
	$line = <SOURCE1>; 

	#ChrC
	$line = <SOURCE1>; 

	#Total
	$line = <SOURCE1>; 
	chomp $line;
	($chrom, $length, $bases, $mean, $min, $max) = split /\s+/, $line; 
	$total_coverage = $mean;


    #----------------------------------------------------------------------
	$source_file2_path = "/zfs/Arabidopsis/work/At_Reseq/bwa_out/".$target_sample."/".$target_sample.".per-base.bed.gz";
	open $fh, "gzip -dc $source_file2_path |";
	$count_base = 0; #for genome coverage X10

	while ($line = <$fh>) {
		chomp $line;
		
		($Chr, $start, $end, $depth) = split /\s+/, $line; 
	
		if($depth > 10){
			$count_base_length = $end - $start; 
			$count_base = $count_base + $count_base_length;
		}
		
	}#close while()

	$genome_cover_rate_x10 = $count_base / $total_base;
	$genome_cover_rate_x10_out = sprintf("%.6f",$genome_cover_rate_x10);
    #----------------------------------------------------------------------


	print OUT $target_sample, "\t", $total_coverage, "\t", $total_nDNA_coverage_out, "\t", $genome_cover_rate_x10_out, "\n";


}

close(OUT);

