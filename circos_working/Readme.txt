#Readme.txt

Step1: Preparation of input files for circos

>perl Make.circos.input.pl < ../../vcf_out/AT.all.list.mutations.txt


Step2: Execution of circos

>sh ./Plot.cirocos.diagram.sh
(circos -conf circos.AT.mutations.conf -outputfile AT.mutations.png)
