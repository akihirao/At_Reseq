#!/usr/bin/perl -i
#ClusterBedBlock.pl
#by HIRAO Akira


#priori setting
$line = <>;
chomp $line;
($pre_Chr, $pre_start, $pre_end, $pre_cluster_ID) = split /\s+/, $line; 
$out_Chr = $pre_Chr;
$out_start = $pre_start; 

$count = 1;

while ($line = <>) {
	chomp $line;

#	print $line, "\n";
	
	($Chr, $start, $end, $cluster_ID) = split /\s+/, $line; 
	
	if($cluster_ID != $pre_cluster_ID){
		$out_end = $pre_end;
		if($count >= 2){
			print $out_Chr, "\t", $out_start, "\t",$out_end, "\n"; 
		}

		$out_Chr = $Chr;
		$out_start = $start;
		$count = 0;

		$pre_Chr = $Chr;
		$pre_start = $start;
		$pre_end = $end;
		$pre_cluster_ID = $cluster_ID; 

	}else{
		$count = $count + 1;

		$pre_Chr = $Chr;
		$pre_start = $start;
		$pre_end = $end;
		$pre_cluster_ID = $cluster_ID; 
	}

		    

}#close while()
