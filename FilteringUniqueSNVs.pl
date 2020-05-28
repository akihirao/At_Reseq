#!/usr/bin/perl -i
#FilteringUniqueSNVs.pl
#from 2019,4.8 to 2019,5.22
#by HIRAO Akira



while ($line = <>) {
	chomp $line;
	if($line =~ m/^##/){
		print $line, "\n";
	}elsif($line =~ m/^#CHROM/) {
		print $line, "\n";

		$line = <>;
		chomp $line;
		($pre2_Chr, $pre2_position, @pre2_info) = split /\s+/, $line; 

		$line = <>;
		chomp $line;
		($pre1_Chr, $pre1_position, @pre1_info) = split /\s+/, $line; 

		$ini_position_dif = $pre1_position - $pre2_position;
	
		if($ini_position_dif > 150){
			print $pre2_Chr, "\t", $pre2_position, "\t";
			$pre2_info_last = pop @pre2_info;
			foreach $output_pre2 (@pre2_info){
				print $output_pre2, "\t";
			}
			print $pre2_info_last, "\n";   	
		}
	}else{
		($Chr, $position, @info) = split /\s+/, $line; 
	    
		$pre_position_dif = abs($pre1_position - $pre2_position);
 		$post_position_dif = abs($position - $pre1_position);

    	if($pre1_Chr ne $pre2_Chr){
    	
    		if($post_position_dif > 150){
	  			
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
    			#no print
    		}

		}else{ #Chr equal

			if($pre_position_dif <= 150 or $post_position_dif <= 150){

				#no print
							
				$pre2_Chr = $pre1_Chr;
				$pre2_position = $pre1_position;
				@pre2_info = @pre1_info;
					
	 			$pre1_Chr = $Chr;
				$pre1_position = $position;
				@pre1_info = @info;

			}else{

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

			}
		}

	}
}#close while()
