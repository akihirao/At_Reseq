#!/usr/bin/perl -i
#Make.circos.input.pl
#by HIRAO Akira
#How to use: perl Make.circos.input.pl < ../../vcf_out/AT.all.list.mutations.txt 

@treat_lab = ("control","low","middle","high");

open (OUT_CON_SBS, ">AT.mutations.$treat_lab[0].sbs.txt");
open (OUT_CON_DEL, ">AT.mutations.$treat_lab[0].deletion.txt");
open (OUT_CON_INS, ">AT.mutations.$treat_lab[0].insertion.txt");

open (OUT_LOW_SBS, ">AT.mutations.$treat_lab[1].sbs.txt");
open (OUT_LOW_DEL, ">AT.mutations.$treat_lab[1].deletion.txt");
open (OUT_LOW_INS, ">AT.mutations.$treat_lab[1].insertion.txt");

open (OUT_MID_SBS, ">AT.mutations.$treat_lab[2].sbs.txt");
open (OUT_MID_DEL, ">AT.mutations.$treat_lab[2].deletion.txt");
open (OUT_MID_INS, ">AT.mutations.$treat_lab[2].insertion.txt");

open (OUT_HIGH_SBS, ">AT.mutations.$treat_lab[3].sbs.txt");
open (OUT_HIGH_DEL, ">AT.mutations.$treat_lab[3].deletion.txt");
open (OUT_HIGH_INS, ">AT.mutations.$treat_lab[3].insertion.txt");

$sbs_color_lab = "color=blue";
$del_color_lab = "color=orange";
$ins_color_lab = "color=green";

$line = <>;

while ($line = <>) {
    chomp $line;
    ($CHROM, $POS, $REF, $ALT, $TYPE, $LENGTH, $SAMPLE1, $ZYGOSITY1, $SAMPLE2, $ZYGOSITY, $SAMPLE3, $ZYGOSITY3, $TREATMENT, $DOSE)  = split /\s+/, $line;
    $CHROM_out = $CHROM;
    substr($CHROM_out,0,1,"c");
    $pre_POS = $POS - 1;

    #Control treatment
    if($TREATMENT eq "Control" && $TYPE eq "SBS"){
        print OUT_CON_SBS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $sbs_color_lab, "\n";
    }
    if($TREATMENT eq "Control" && $TYPE eq "Deletion"){
        print OUT_CON_DEL $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $del_color_lab, "\n";
    }
    if($TREATMENT eq "Control" && $TYPE eq "Insertion"){
        print OUT_CON_INS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $ins_color_lab, "\n";
    }

    #Low treatment
    if($TREATMENT eq "Low" && $TYPE eq "SBS"){
        print OUT_LOW_SBS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $sbs_color_lab, "\n";
    }
    if($TREATMENT eq "Low" && $TYPE eq "Deletion"){
        print OUT_LOW_DEL $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $del_color_lab, "\n";
    }
    if($TREATMENT eq "Low" && $TYPE eq "Insertion"){
        print OUT_LOW_INS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $ins_color_lab, "\n";
    }

    #Middle treatment
    if($TREATMENT eq "Middle" && $TYPE eq "SBS"){
        print OUT_MID_SBS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $sbs_color_lab, "\n";
    }
    if($TREATMENT eq "Middle" && $TYPE eq "Deletion"){
        print OUT_MID_DEL $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $del_color_lab, "\n";
    }
    if($TREATMENT eq "Middle" && $TYPE eq "Insertion"){
        print OUT_MID_INS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $ins_color_lab, "\n";
    }

    #High treatment
    if($TREATMENT eq "High" && $TYPE eq "SBS"){
        print OUT_HIGH_SBS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $sbs_color_lab, "\n";
    }
    if($TREATMENT eq "High" && $TYPE eq "Deletion"){
        print OUT_HIGH_DEL $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $del_color_lab, "\n";
    }
    if($TREATMENT eq "High" && $TYPE eq "Insertion"){
        print OUT_HIGH_INS $CHROM_out, "\t", $pre_POS, "\t", $POS, "\t", "+", "\t", $ins_color_lab, "\n";
    }


}

close (OUT_CON_SBS);
close (OUT_CON_DEL);
close (OUT_CON_INS);

close (OUT_LOW_SBS);
close (OUT_LOW_DEL);
close (OUT_LOW_INS);

close (OUT_MID_SBS);
close (OUT_MID_DEL);
close (OUT_MID_INS);

close (OUT_HIGH_SBS);
close (OUT_HIGH_DEL);
close (OUT_HIGH_INS);
