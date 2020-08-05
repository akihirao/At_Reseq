#snpEff_effect_summary.R
#by Akira Hirao


args <- commandArgs()
sample.ID <- args[6]

system("cat snpEff_genes.txt | tail -n +2 > snpEff_genes.edited.txt")
system("sed -i -e '1,1s:^#::' snpEff_genes.edited.txt")

snpEff_genes.data <- read.table("snpEff_genes.edited.txt", header=T, sep="\t")
high.effect <- sum(snpEff_genes.data$variants_impact_HIGH)
moderate.effect <- sum(snpEff_genes.data$variants_impact_MODERATE)
modifier.effect <- sum(snpEff_genes.data$variants_impact_MODIFIER)
low.effect <- sum(snpEff_genes.data$variants_impact_LOW)

output.file.name <- paste(sample.ID, "snpEff_effect_summary.txt",sep=".")
output.contents <- data.frame(Sample=sample.ID, High=high.effect, Moderate=moderate.effect, Modifier=modifier.effect, Low=low.effect)
#writeLines(output.header, output.file.name, sep="\t")
write.table(output.contents, output.file.name, row.names=F,quote=F)