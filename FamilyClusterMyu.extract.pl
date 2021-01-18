#!/usr/bin/perl -i
#FamilyClusterMyu.extract.pl
#by HIRAO Akira


@family_lab = ("A010","A020","A030","A110","A120","A130","A210","A220","A230","A310","A320","A330");

open (OUT_FAMILY_ALL, ">AT48.family.clustered.mu.vcf");

open (OUT_A010, ">$family_lab[0].family.vcf");
open (OUT_A020, ">$family_lab[1].family.vcf");
open (OUT_A030, ">$family_lab[2].family.vcf");

open (OUT_A110, ">$family_lab[3].family.vcf");
open (OUT_A120, ">$family_lab[4].family.vcf");
open (OUT_A130, ">$family_lab[5].family.vcf");

open (OUT_A210, ">$family_lab[6].family.vcf");
open (OUT_A220, ">$family_lab[7].family.vcf");
open (OUT_A230, ">$family_lab[8].family.vcf");

open (OUT_A310, ">$family_lab[9].family.vcf");
open (OUT_A320, ">$family_lab[10].family.vcf");
open (OUT_A330, ">$family_lab[11].family.vcf");

open (OUT_A010_ALL, ">$family_lab[0].family.all.vcf");
open (OUT_A020_ALL, ">$family_lab[1].family.all.vcf");
open (OUT_A030_ALL, ">$family_lab[2].family.all.vcf");

open (OUT_A110_ALL, ">$family_lab[3].family.all.vcf");
open (OUT_A120_ALL, ">$family_lab[4].family.all.vcf");
open (OUT_A130_ALL, ">$family_lab[5].family.all.vcf");

open (OUT_A210_ALL, ">$family_lab[6].family.all.vcf");
open (OUT_A220_ALL, ">$family_lab[7].family.all.vcf");
open (OUT_A230_ALL, ">$family_lab[8].family.all.vcf");

open (OUT_A310_ALL, ">$family_lab[9].family.all.vcf");
open (OUT_A320_ALL, ">$family_lab[10].family.all.vcf");
open (OUT_A330_ALL, ">$family_lab[11].family.all.vcf");


while ($line = <>) {
    chomp $line;
    if($line !~ m/^#/){        
        #($Chr, $end, @info) = split /\s+/, $line;
        ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $A011, $A012, $A014, $A021, $A023, $A024, $A031, $A033, $A034, $A113, $A114, $A115, $A121, $A123, $A125, $A131, $A132, $A135, $A212, $A213, $A215, $A221, $A222, $A225, $A231, $A233, $A234, $A311, $A313, $A314, $A321, $A323, $A324, $A331, $A332, $A334, $AT01, $AT02, $AT03, $AT11, $AT12, $AT13, $AT21, $AT22, $AT23, $AT31, $AT32, $AT33)  = split /\s+/, $line;

        ($GT_A011, @info) = split /:/, $A011;
        ($GT_A012, @info) = split /:/, $A012;
        ($GT_A014, @info) = split /:/, $A014;
        ($GT_A021, @info) = split /:/, $A021;
        ($GT_A023, @info) = split /:/, $A023;
        ($GT_A024, @info) = split /:/, $A024;
        ($GT_A031, @info) = split /:/, $A031;
        ($GT_A033, @info) = split /:/, $A033;
        ($GT_A034, @info) = split /:/, $A034;
        ($GT_A113, @info) = split /:/, $A113;
        ($GT_A114, @info) = split /:/, $A114;
        ($GT_A115, @info) = split /:/, $A115;
        ($GT_A121, @info) = split /:/, $A121;
        ($GT_A123, @info) = split /:/, $A123;
        ($GT_A125, @info) = split /:/, $A125;
        ($GT_A131, @info) = split /:/, $A131;
        ($GT_A132, @info) = split /:/, $A132;
        ($GT_A135, @info) = split /:/, $A135;
        ($GT_A212, @info) = split /:/, $A212;
        ($GT_A213, @info) = split /:/, $A213;
        ($GT_A215, @info) = split /:/, $A215;
        ($GT_A221, @info) = split /:/, $A221;
        ($GT_A222, @info) = split /:/, $A222;
        ($GT_A225, @info) = split /:/, $A225;
        ($GT_A231, @info) = split /:/, $A231;
        ($GT_A233, @info) = split /:/, $A233;
        ($GT_A234, @info) = split /:/, $A234;
        ($GT_A311, @info) = split /:/, $A311;
        ($GT_A313, @info) = split /:/, $A313;
        ($GT_A314, @info) = split /:/, $A314;
        ($GT_A321, @info) = split /:/, $A321;
        ($GT_A323, @info) = split /:/, $A323;
        ($GT_A324, @info) = split /:/, $A324;
        ($GT_A331, @info) = split /:/, $A331;
        ($GT_A332, @info) = split /:/, $A332;
        ($GT_A334, @info) = split /:/, $A334;
        ($GT_AT01, @info) = split /:/, $AT01;
        ($GT_AT02, @info) = split /:/, $AT02;
        ($GT_AT03, @info) = split /:/, $AT03;
        ($GT_AT11, @info) = split /:/, $AT11;
        ($GT_AT12, @info) = split /:/, $AT12;
        ($GT_AT13, @info) = split /:/, $AT13;
        ($GT_AT21, @info) = split /:/, $AT21;
        ($GT_AT22, @info) = split /:/, $AT22;
        ($GT_AT23, @info) = split /:/, $AT23;
        ($GT_AT31, @info) = split /:/, $AT31;
        ($GT_AT32, @info) = split /:/, $AT32;
        ($GT_AT33, @info) = split /:/, $AT33;

        $GT_A011_a = substr($GT_A011,0,1); $GT_A011_b = substr($GT_A011,2,1); $GT_sum_A011 = $GT_A011_a + $GT_A011_b;if($GT_A011_b > 0) {$GT_A011_Alt = 1}else{$GT_A011_Alt = 0};
        $GT_A012_a = substr($GT_A012,0,1); $GT_A012_b = substr($GT_A012,2,1); $GT_sum_A012 = $GT_A012_a + $GT_A012_b;if($GT_A012_b > 0) {$GT_A012_Alt = 1}else{$GT_A012_Alt = 0};
        $GT_A014_a = substr($GT_A014,0,1); $GT_A014_b = substr($GT_A014,2,1); $GT_sum_A014 = $GT_A014_a + $GT_A014_b;if($GT_A014_b > 0) {$GT_A014_Alt = 1}else{$GT_A014_Alt = 0};
        $GT_AT01_a = substr($GT_AT01,0,1); $GT_AT01_b = substr($GT_AT01,2,1); $GT_sum_AT01 = $GT_AT01_a + $GT_AT01_b;if($GT_AT01_b > 0) {$GT_AT01_Alt = 1}else{$GT_AT01_Alt = 0};

        $GT_A021_a = substr($GT_A021,0,1); $GT_A021_b = substr($GT_A021,2,1); $GT_sum_A021 = $GT_A021_a + $GT_A021_b;if($GT_A021_b > 0) {$GT_A021_Alt = 1}else{$GT_A021_Alt = 0};
        $GT_A023_a = substr($GT_A023,0,1); $GT_A023_b = substr($GT_A023,2,1); $GT_sum_A023 = $GT_A023_a + $GT_A023_b;if($GT_A023_b > 0) {$GT_A023_Alt = 1}else{$GT_A023_Alt = 0};
        $GT_A024_a = substr($GT_A024,0,1); $GT_A024_b = substr($GT_A024,2,1); $GT_sum_A024 = $GT_A024_a + $GT_A024_b;if($GT_A024_b > 0) {$GT_A024_Alt = 1}else{$GT_A024_Alt = 0};
        $GT_AT02_a = substr($GT_AT02,0,1); $GT_AT02_b = substr($GT_AT02,2,1); $GT_sum_AT02 = $GT_AT02_a + $GT_AT02_b;if($GT_AT02_b > 0) {$GT_AT02_Alt = 1}else{$GT_AT02_Alt = 0};

        $GT_A031_a = substr($GT_A031,0,1); $GT_A031_b = substr($GT_A031,2,1); $GT_sum_A031 = $GT_A031_a + $GT_A031_b;if($GT_A031_b > 0) {$GT_A031_Alt = 1}else{$GT_A031_Alt = 0};
        $GT_A033_a = substr($GT_A033,0,1); $GT_A033_b = substr($GT_A033,2,1); $GT_sum_A033 = $GT_A033_a + $GT_A033_b;if($GT_A033_b > 0) {$GT_A033_Alt = 1}else{$GT_A033_Alt = 0};
        $GT_A034_a = substr($GT_A034,0,1); $GT_A034_b = substr($GT_A034,2,1); $GT_sum_A034 = $GT_A034_a + $GT_A034_b;if($GT_A034_b > 0) {$GT_A034_Alt = 1}else{$GT_A034_Alt = 0};
        $GT_AT03_a = substr($GT_AT03,0,1); $GT_AT03_b = substr($GT_AT03,2,1); $GT_sum_AT03 = $GT_AT03_a + $GT_AT03_b;if($GT_AT03_b > 0) {$GT_AT03_Alt = 1}else{$GT_AT03_Alt = 0};

        $GT_A113_a = substr($GT_A113,0,1); $GT_A113_b = substr($GT_A113,2,1); $GT_sum_A113 = $GT_A113_a + $GT_A113_b;if($GT_A113_b > 0) {$GT_A113_Alt = 1}else{$GT_A113_Alt = 0};
        $GT_A114_a = substr($GT_A114,0,1); $GT_A114_b = substr($GT_A114,2,1); $GT_sum_A114 = $GT_A114_a + $GT_A114_b;if($GT_A114_b > 0) {$GT_A114_Alt = 1}else{$GT_A114_Alt = 0};
        $GT_A115_a = substr($GT_A115,0,1); $GT_A115_b = substr($GT_A115,2,1); $GT_sum_A115 = $GT_A115_a + $GT_A115_b;if($GT_A115_b > 0) {$GT_A115_Alt = 1}else{$GT_A115_Alt = 0};
        $GT_AT11_a = substr($GT_AT11,0,1); $GT_AT11_b = substr($GT_AT11,2,1); $GT_sum_AT11 = $GT_AT11_a + $GT_AT11_b;if($GT_AT11_b > 0) {$GT_AT11_Alt = 1}else{$GT_AT11_Alt = 0};

        $GT_A121_a = substr($GT_A121,0,1); $GT_A121_b = substr($GT_A121,2,1); $GT_sum_A121 = $GT_A121_a + $GT_A121_b;if($GT_A121_b > 0) {$GT_A121_Alt = 1}else{$GT_A121_Alt = 0};
        $GT_A123_a = substr($GT_A123,0,1); $GT_A123_b = substr($GT_A123,2,1); $GT_sum_A123 = $GT_A123_a + $GT_A123_b;if($GT_A123_b > 0) {$GT_A123_Alt = 1}else{$GT_A123_Alt = 0};
        $GT_A125_a = substr($GT_A125,0,1); $GT_A125_b = substr($GT_A125,2,1); $GT_sum_A125 = $GT_A125_a + $GT_A125_b;if($GT_A125_b > 0) {$GT_A125_Alt = 1}else{$GT_A125_Alt = 0};
        $GT_AT12_a = substr($GT_AT12,0,1); $GT_AT12_b = substr($GT_AT12,2,1); $GT_sum_AT12 = $GT_AT12_a + $GT_AT12_b;if($GT_AT12_b > 0) {$GT_AT12_Alt = 1}else{$GT_AT12_Alt = 0};

        $GT_A131_a = substr($GT_A131,0,1); $GT_A131_b = substr($GT_A131,2,1); $GT_sum_A131 = $GT_A131_a + $GT_A131_b;if($GT_A131_b > 0) {$GT_A131_Alt = 1}else{$GT_A131_Alt = 0};
        $GT_A132_a = substr($GT_A132,0,1); $GT_A132_b = substr($GT_A132,2,1); $GT_sum_A132 = $GT_A132_a + $GT_A132_b;if($GT_A132_b > 0) {$GT_A132_Alt = 1}else{$GT_A132_Alt = 0};
        $GT_A135_a = substr($GT_A135,0,1); $GT_A135_b = substr($GT_A135,2,1); $GT_sum_A135 = $GT_A135_a + $GT_A135_b;if($GT_A135_b > 0) {$GT_A135_Alt = 1}else{$GT_A135_Alt = 0};
        $GT_AT13_a = substr($GT_AT13,0,1); $GT_AT13_b = substr($GT_AT13,2,1); $GT_sum_AT13 = $GT_AT13_a + $GT_AT13_b;if($GT_AT13_b > 0) {$GT_AT13_Alt = 1}else{$GT_AT13_Alt = 0};

        $GT_A212_a = substr($GT_A212,0,1); $GT_A212_b = substr($GT_A212,2,1); $GT_sum_A212 = $GT_A212_a + $GT_A212_b;if($GT_A212_b > 0) {$GT_A212_Alt = 1}else{$GT_A212_Alt = 0};
        $GT_A213_a = substr($GT_A213,0,1); $GT_A213_b = substr($GT_A213,2,1); $GT_sum_A213 = $GT_A213_a + $GT_A213_b;if($GT_A213_b > 0) {$GT_A213_Alt = 1}else{$GT_A213_Alt = 0};
        $GT_A215_a = substr($GT_A215,0,1); $GT_A215_b = substr($GT_A215,2,1); $GT_sum_A215 = $GT_A215_a + $GT_A215_b;if($GT_A215_b > 0) {$GT_A215_Alt = 1}else{$GT_A215_Alt = 0};
        $GT_AT21_a = substr($GT_AT21,0,1); $GT_AT21_b = substr($GT_AT21,2,1); $GT_sum_AT21 = $GT_AT21_a + $GT_AT21_b;if($GT_AT21_b > 0) {$GT_AT21_Alt = 1}else{$GT_AT21_Alt = 0};

        $GT_A221_a = substr($GT_A221,0,1); $GT_A221_b = substr($GT_A221,2,1); $GT_sum_A221 = $GT_A221_a + $GT_A221_b;if($GT_A221_b > 0) {$GT_A221_Alt = 1}else{$GT_A221_Alt = 0};
        $GT_A222_a = substr($GT_A222,0,1); $GT_A222_b = substr($GT_A222,2,1); $GT_sum_A222 = $GT_A222_a + $GT_A222_b;if($GT_A222_b > 0) {$GT_A222_Alt = 1}else{$GT_A222_Alt = 0};
        $GT_A225_a = substr($GT_A225,0,1); $GT_A225_b = substr($GT_A225,2,1); $GT_sum_A225 = $GT_A225_a + $GT_A225_b;if($GT_A225_b > 0) {$GT_A225_Alt = 1}else{$GT_A225_Alt = 0};
        $GT_AT22_a = substr($GT_AT22,0,1); $GT_AT22_b = substr($GT_AT22,2,1); $GT_sum_AT22 = $GT_AT22_a + $GT_AT22_b;if($GT_AT22_b > 0) {$GT_AT22_Alt = 1}else{$GT_AT22_Alt = 0};

        $GT_A231_a = substr($GT_A231,0,1); $GT_A231_b = substr($GT_A231,2,1); $GT_sum_A231 = $GT_A231_a + $GT_A231_b;if($GT_A231_b > 0) {$GT_A231_Alt = 1}else{$GT_A231_Alt = 0};
        $GT_A233_a = substr($GT_A233,0,1); $GT_A233_b = substr($GT_A233,2,1); $GT_sum_A233 = $GT_A233_a + $GT_A233_b;if($GT_A233_b > 0) {$GT_A233_Alt = 1}else{$GT_A233_Alt = 0};
        $GT_A234_a = substr($GT_A234,0,1); $GT_A234_b = substr($GT_A234,2,1); $GT_sum_A234 = $GT_A234_a + $GT_A234_b;if($GT_A234_b > 0) {$GT_A234_Alt = 1}else{$GT_A234_Alt = 0};
        $GT_AT23_a = substr($GT_AT23,0,1); $GT_AT23_b = substr($GT_AT23,2,1); $GT_sum_AT23 = $GT_AT23_a + $GT_AT23_b;if($GT_AT23_b > 0) {$GT_AT23_Alt = 1}else{$GT_AT23_Alt = 0};

        $GT_A311_a = substr($GT_A311,0,1); $GT_A311_b = substr($GT_A311,2,1); $GT_sum_A311 = $GT_A311_a + $GT_A311_b;if($GT_A311_b > 0) {$GT_A311_Alt = 1}else{$GT_A311_Alt = 0};
        $GT_A313_a = substr($GT_A313,0,1); $GT_A313_b = substr($GT_A313,2,1); $GT_sum_A313 = $GT_A313_a + $GT_A313_b;if($GT_A313_b > 0) {$GT_A313_Alt = 1}else{$GT_A313_Alt = 0};
        $GT_A314_a = substr($GT_A314,0,1); $GT_A314_b = substr($GT_A314,2,1); $GT_sum_A314 = $GT_A314_a + $GT_A314_b;if($GT_A314_b > 0) {$GT_A314_Alt = 1}else{$GT_A314_Alt = 0};
        $GT_AT31_a = substr($GT_AT31,0,1); $GT_AT31_b = substr($GT_AT31,2,1); $GT_sum_AT31 = $GT_AT31_a + $GT_AT31_b;if($GT_AT31_b > 0) {$GT_AT31_Alt = 1}else{$GT_AT31_Alt = 0};

        $GT_A321_a = substr($GT_A321,0,1); $GT_A321_b = substr($GT_A321,2,1); $GT_sum_A321 = $GT_A321_a + $GT_A321_b;if($GT_A321_b > 0) {$GT_A321_Alt = 1}else{$GT_A321_Alt = 0};
        $GT_A323_a = substr($GT_A323,0,1); $GT_A323_b = substr($GT_A323,2,1); $GT_sum_A323 = $GT_A323_a + $GT_A323_b;if($GT_A323_b > 0) {$GT_A323_Alt = 1}else{$GT_A323_Alt = 0};
        $GT_A324_a = substr($GT_A324,0,1); $GT_A324_b = substr($GT_A324,2,1); $GT_sum_A324 = $GT_A324_a + $GT_A324_b;if($GT_A324_b > 0) {$GT_A324_Alt = 1}else{$GT_A324_Alt = 0};
        $GT_AT32_a = substr($GT_AT32,0,1); $GT_AT32_b = substr($GT_AT32,2,1); $GT_sum_AT32 = $GT_AT32_a + $GT_AT32_b;if($GT_AT32_b > 0) {$GT_AT32_Alt = 1}else{$GT_AT32_Alt = 0};

        $GT_A331_a = substr($GT_A331,0,1); $GT_A331_b = substr($GT_A331,2,1); $GT_sum_A331 = $GT_A331_a + $GT_A331_b;if($GT_A331_b > 0) {$GT_A331_Alt = 1}else{$GT_A331_Alt = 0};
        $GT_A332_a = substr($GT_A332,0,1); $GT_A332_b = substr($GT_A332,2,1); $GT_sum_A332 = $GT_A332_a + $GT_A332_b;if($GT_A332_b > 0) {$GT_A332_Alt = 1}else{$GT_A332_Alt = 0};
        $GT_A334_a = substr($GT_A334,0,1); $GT_A334_b = substr($GT_A334,2,1); $GT_sum_A334 = $GT_A334_a + $GT_A334_b;if($GT_A334_b > 0) {$GT_A334_Alt = 1}else{$GT_A334_Alt = 0};
        $GT_AT33_a = substr($GT_AT33,0,1); $GT_AT33_b = substr($GT_AT33,2,1); $GT_sum_AT33 = $GT_AT33_a + $GT_AT33_b;if($GT_AT33_b > 0) {$GT_AT33_Alt = 1}else{$GT_AT33_Alt = 0};

#        print $GT_A011_Alt, "\n";

        @GT_48indiv_vec = ($GT_sum_A011, $GT_sum_A012, $GT_sum_A014, $GT_sum_A021, $GT_sum_A023, $GT_sum_A024, $GT_sum_A031, $GT_sum_A033, $GT_sum_A034, $GT_sum_A113, $GT_sum_A114, $GT_sum_A115, $GT_sum_A121, $GT_sum_A123, $GT_sum_A125, $GT_sum_A131, $GT_sum_A132, $GT_sum_A135, $GT_sum_A212, $GT_sum_A213, $GT_sum_A215, $GT_sum_A221, $GT_sum_A222, $GT_sum_A225, $GT_sum_A231, $GT_sum_A233, $GT_sum_A234, $GT_sum_A311, $GT_sum_A313, $GT_sum_A314, $GT_sum_A321, $GT_sum_A323, $GT_sum_A324, $GT_sum_A331, $GT_sum_A332, $GT_sum_A334, $GT_sum_AT01, $GT_sum_AT02, $GT_sum_AT03, $GT_sum_AT11, $GT_sum_AT12, $GT_sum_AT13, $GT_sum_AT21, $GT_sum_AT22, $GT_sum_AT23, $GT_sum_AT31, $GT_sum_AT32, $GT_sum_AT33);
        @GT_48indiv_alt = ($GT_A011_Alt, $GT_A012_Alt, $GT_A014_Alt, $GT_A021_Alt, $GT_A023_Alt, $GT_A024_Alt, $GT_A031_Alt, $GT_A033_Alt, $GT_A034_Alt, $GT_A113_Alt, $GT_A114_Alt, $GT_A115_Alt, $GT_A121_Alt, $GT_A123_Alt, $GT_A125_Alt, $GT_A131_Alt, $GT_A132_Alt, $GT_A135_Alt, $GT_A212_Alt, $GT_A213_Alt, $GT_A215_Alt, $GT_A221_Alt, $GT_A222_Alt, $GT_A225_Alt, $GT_A231_Alt, $GT_A233_Alt, $GT_A234_Alt, $GT_A311_Alt, $GT_A313_Alt, $GT_A314_Alt, $GT_A321_Alt, $GT_A323_Alt, $GT_A324_Alt, $GT_A331_Alt, $GT_A332_Alt, $GT_A334_Alt, $GT_AT01_Alt, $GT_AT02_Alt, $GT_AT03_Alt, $GT_AT11_Alt, $GT_AT12_Alt, $GT_AT13_Alt, $GT_AT21_Alt, $GT_AT22_Alt, $GT_AT23_Alt, $GT_AT31_Alt, $GT_AT32_Alt, $GT_AT33_Alt);
    

        $A010_mu = 0; $non_A010_mu = 0;  #initializing zero
        @A010_non = @GT_48indiv_alt; @A010_fam = splice (@A010_non,0,3);$A010_mu += $_ for @A010_fam;$non_A010_mu += $_ for @A010_non;
        if($A010_mu >= 2 and $non_A010_mu == 0 and $AT010 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A010 $line, "\n";

        }
         if($A010_mu == 3 and $non_A010_mu == 0 and $AT010 == 0){
            print OUT_A010_ALL $line, "\n";
        }
  
  
        $A020_mu = 0; $non_A020_mu = 0;  #initializing zero
        @A020_non = @GT_48indiv_alt; @A020_fam = splice (@A020_non,3,3);$A020_mu += $_ for @A020_fam;$non_A020_mu += $_ for @A020_non;
        if($A020_mu >= 2 and $non_A020_mu == 0 and $AT020 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A020 $line, "\n";
        }
         if($A020_mu == 3 and $non_A020_mu == 0 and $AT020 == 0){
            print OUT_A020_ALL $line, "\n";
        }
 
        $A030_mu = 0; $non_A030_mu = 0;  #initializing zero
        @A030_non = @GT_48indiv_alt; @A030_fam = splice (@A030_non,6,3);$A030_mu += $_ for @A030_fam;$non_A030_mu += $_ for @A030_non;
        if($A030_mu >= 2 and $non_A030_mu == 0 and $AT030 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A030 $line, "\n";
        }
         if($A030_mu == 3 and $non_A030_mu == 0 and $AT030 == 0){
            print OUT_A030_ALL $line, "\n";
        }
 
        $A110_mu = 0; $non_A110_mu = 0;  #initializing zero
        @A110_non = @GT_48indiv_alt; @A110_fam = splice (@A110_non,9,3);$A110_mu += $_ for @A110_fam;$non_A110_mu += $_ for @A110_non;
        if($A110_mu >= 2 and $non_A110_mu == 0 and $AT110 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A110 $line, "\n";
        }
         if($A110_mu == 3 and $non_A110_mu == 0 and $AT110 == 0){
            print OUT_A110_ALL $line, "\n";
        }
 
        $A120_mu = 0; $non_A120_mu = 0;  #initializing zero
        @A120_non = @GT_48indiv_alt; @A120_fam = splice (@A120_non,12,3);$A120_mu += $_ for @A120_fam;$non_A120_mu += $_ for @A120_non;
        if($A120_mu >= 2 and $non_A120_mu == 0 and $AT120 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A120 $line, "\n";
        }
         if($A120_mu == 3 and $non_A120_mu == 0 and $AT120 == 0){
            print OUT_A120_ALL $line, "\n";
        }
 
        $A130_mu = 0; $non_A130_mu = 0;  #initializing zero
        @A130_non = @GT_48indiv_alt; @A130_fam = splice (@A130_non,15,3);$A130_mu += $_ for @A130_fam;$non_A130_mu += $_ for @A130_non;
        if($A130_mu >= 2 and $non_A130_mu == 0 and $AT130 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A130 $line, "\n";
        }
         if($A130_mu == 3 and $non_A130_mu == 0 and $AT130 == 0){
            print OUT_A130_ALL $line, "\n";
        }
 
        $A210_mu = 0; $non_A210_mu = 0;  #initializing zero
        @A210_non = @GT_48indiv_alt; @A210_fam = splice (@A210_non,18,3);$A210_mu += $_ for @A210_fam;$non_A210_mu += $_ for @A210_non;
        if($A210_mu >= 2 and $non_A210_mu == 0 and $AT210 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A210 $line, "\n";
        }
         if($A210_mu == 3 and $non_A210_mu == 0 and $AT210 == 0){
            print OUT_A210_ALL $line, "\n";
        }
 
        $A220_mu = 0; $non_A220_mu = 0;  #initializing zero
        @A220_non = @GT_48indiv_alt; @A220_fam = splice (@A220_non,21,3);$A220_mu += $_ for @A220_fam;$non_A220_mu += $_ for @A220_non;
        if($A220_mu >= 2 and $non_A220_mu == 0 and $AT220 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A220 $line, "\n";
        }
         if($A220_mu == 3 and $non_A220_mu == 0 and $AT220 == 0){
            print OUT_A220_ALL $line, "\n";
        }
 
        $A230_mu = 0; $non_A230_mu = 0;  #initializing zero
        @A230_non = @GT_48indiv_alt; @A230_fam = splice (@A230_non,24,3);$A230_mu += $_ for @A230_fam;$non_A230_mu += $_ for @A230_non;
        if($A230_mu >= 2 and $non_A230_mu == 0 and $AT230 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A230 $line, "\n";
        }
         if($A230_mu == 3 and $non_A230_mu == 0 and $AT230 == 0){
            print OUT_A230_ALL $line, "\n";
        }
 
        $A310_mu = 0; $non_A310_mu = 0;  #initializing zero
        @A310_non = @GT_48indiv_alt; @A310_fam = splice (@A310_non,27,3);$A310_mu += $_ for @A310_fam;$non_A310_mu += $_ for @A310_non;
        if($A310_mu >= 2 and $non_A310_mu == 0 and $AT310 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A310 $line, "\n";
        }       
         if($A310_mu == 3 and $non_A310_mu == 0 and $AT310 == 0){
            print OUT_A310_ALL $line, "\n";
        }
 
        $A320_mu = 0; $non_A320_mu = 0;  #initializing zero
        @A320_non = @GT_48indiv_alt; @A320_fam = splice (@A320_non,30,3);$A320_mu += $_ for @A320_fam;$non_A320_mu += $_ for @A320_non;
        if($A320_mu >= 2 and $non_A320_mu == 0 and $AT320 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A320 $line, "\n";
        }
          if($A320_mu == 3 and $non_A320_mu == 0 and $AT320 == 0){
            print OUT_A320_ALL $line, "\n";
        }
 
        $A330_mu = 0; $non_A330_mu = 0;  #initializing zero
        @A330_non = @GT_48indiv_alt; @A330_fam = splice (@A330_non,33,3);$A330_mu += $_ for @A330_fam;$non_A330_mu += $_ for @A330_non;
        if($A330_mu >= 2 and $non_A330_mu == 0 and $AT330 == 0){
            print OUT_FAMILY_ALL $line, "\n";
            print OUT_A330 $line, "\n";
        }
         if($A330_mu == 3 and $non_A330_mu == 0 and $AT330 == 0){
            print OUT_A330_ALL $line, "\n";
        }
 
    }else{

        print OUT_FAMILY_ALL $line, "\n";

        print OUT_A010 $line, "\n";
        print OUT_A020 $line, "\n";
        print OUT_A030 $line, "\n";

        print OUT_A110 $line, "\n";
        print OUT_A120 $line, "\n";
        print OUT_A130 $line, "\n";

        print OUT_A210 $line, "\n";
        print OUT_A220 $line, "\n";
        print OUT_A230 $line, "\n";

        print OUT_A310 $line, "\n";
        print OUT_A320 $line, "\n";
        print OUT_A330 $line, "\n";

        print OUT_A010_ALL $line, "\n";
        print OUT_A020_ALL $line, "\n";
        print OUT_A030_ALL $line, "\n";

        print OUT_A110_ALL $line, "\n";
        print OUT_A120_ALL $line, "\n";
        print OUT_A130_ALL $line, "\n";

        print OUT_A210_ALL $line, "\n";
        print OUT_A220_ALL $line, "\n";
        print OUT_A230_ALL $line, "\n";

        print OUT_A310_ALL $line, "\n";
        print OUT_A320_ALL $line, "\n";
        print OUT_A330_ALL $line, "\n";

    }

}

close (OUT_A010);
close (OUT_A020);
close (OUT_A030);

close (OUT_A110);
close (OUT_A120);
close (OUT_A130);

close (OUT_A210);
close (OUT_A220);
close (OUT_A230);

close (OUT_A310);
close (OUT_A320);
close (OUT_A330);
