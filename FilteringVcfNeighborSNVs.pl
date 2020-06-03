#!/usr/bin/perl -i
#FilteringVcfNeighborSNVs.pl
#by HIRAO Akira

$neighbor_threshold = 300; #

#priori setting
$pre2_Chr = "Chr1";
$pre1_Chr = "Chr1";
$pre2_position = 0; 
$pre1_position = 1; 

while ($line = <>) {
	chomp $line;

	if(eof()){
		($last_Chr, $last_position, @last_info) = split /\s+/, $line; 
	}elsif($line =~ m/^#/){
		print $line, "\n";
	}else{
		($Chr, $position, @info) = split /\s+/, $line; 
	    
		$pre_position_dif = abs($pre1_position - $pre2_position);
 		$post_position_dif = abs($position - $pre1_position);


		if($pre_position_dif > $neighbor_threshold and $post_position_dif > $neighbor_threshold){

		 	print $pre1_Chr, "\t", $pre1_position, "\t";

			$pre1_info_last = pop @pre1_info;
			foreach $output_pre1 (@pre1_info){
				print $output_pre1, "\t";
			}
			print $pre1_info_last, "\n";    	

			$pre2_Chr = $pre1_Chr;
			$pre2_position = $pre1_position;
			@pre2_info = @pre1_info;
					
	 		$pre1_Chr = $Chr;
			$pre1_position = $position;
			@pre1_info = @info;

		}else{

			$pre2_Chr = $pre1_Chr;
			$pre2_position = $pre1_position;
			@pre2_info = @pre1_info;
					
	 		$pre1_Chr = $Chr;
			$pre1_position = $position;
			@pre1_info = @info;
			
		}				

	}

}#close while()


#handling final two lines
$last_pre_position_dif = abs($position - $pre2_position); # pre1_position had replaced with pre2_position
$last_post_position_dif = abs($last_position - $position);

if($last_pre_position_dif > $neighbor_threshold and $last_post_position_dif > $neighbor_threshold){
	print $Chr, "\t", $position, "\t";
	$info_last = pop @info;
	foreach $output_info (@info){
		print $output_info, "\t";
	}
	print $info_last, "\n";    	
}

if($last_post_position_dif > $neighbor_threshold){
	print $last_Chr, "\t", $last_position, "\t";
	$last_info_last = pop @last_info;
	foreach $output_last_info (@last_info){
		print $output_last_info, "\t";
	}
	print $last_info_last, "\n";    	
}
