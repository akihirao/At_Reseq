#summarizing_M1M2All.unique_snps.r
#2019,7.19
args <- commandArgs() #vector of parameters
category <- args[3]
sample1 <- args[4]
sample2 <- args[5]
sample3 <- args[6]
sample4 <- args[7]
sample5 <- args[8]
sample6 <- args[9]
sample7 <- args[10]
sample8 <- args[11]
sample9 <- args[12]
sample10 <- args[13]
sample11 <- args[14]
sample12 <- args[15]


input_snps_file <- paste(category,".family.unique_het_snps_list.txt",sep="")
family.snps.list <- read.table(input_snps_file)
no.total.snp <- length(family.snps.list[,1])

family.snps.list$V4 <- factor(family.snps.list$V4, 
	levels=c(category, sample1, sample2, sample3))
family.no_unique_snps <- table(family.snps.list$V4)

#print(family.no_unique_snps)

if(FALSE){
input_indels_file <- paste(category,".family.unique_het_indels_list.txt",sep="")
family.indels.list <- read.table(input_indels_file)
no.total.indel <- length(family.indels.list[,1])

family.indels.list$V4 <- factor(family.indels.list$V4, 
	levels=c(category, sample1, sample2, sample3))
family.no_unique_indels <- table(family.indels.list$V4)

#print(family.no_unique_indels)
}

input_snps_file.filtered <- paste(category,".family.unique_het_snps_list.filtered.txt",sep="")
family.snps.list.filtered <- read.table(input_snps_file.filtered)
no.total.snp.filtered <- length(family.snps.list.filtered[,1])

family.snps.list.filtered$V4 <- factor(family.snps.list.filtered$V4, 
	levels=c(category, sample1, sample2, sample3))
family.no_unique_snps.filtered <- table(family.snps.list.filtered$V4)

#print(family.no_unique_snps.filtered)

if(FALSE){
input_indels_file.filtered <- paste(category,".family.unique_het_indels_list.filtered.txt",sep="")
family.indels.list.filtered <- read.table(input_indels_file.filtered)
no.total.indel.filtered <- length(family.indels.list.filtered[,1])

family.indels.list.filtered$V4 <- factor(family.indels.list.filtered$V4, 
	levels=c(category, sample1, sample2, sample3))
family.no_unique_indels.filtered <- table(family.indels.list.filtered$V4)

#print(family.no_unique_indels.filtered)
}

#output.no_unique.snv <- rbind(family.no_unique_snps, family.no_unique_indels,family.no_unique_snps.filtered,family.no_unique_indels.filtered)
#row.names(output.no_unique.snv) <- c("no.unique.snps","no.unique.indels","no.unique.snps.filtered","no.unique.indels.filtered")

output.no_unique.snv <- rbind(family.no_unique_snps, family.no_unique_snps.filtered)
row.names(output.no_unique.snv) <- c("no.unique.snps","no.unique.snps.filtered")


print(output.no_unique.snv)
write.table(output.no_unique.snv, "summary.no_snv.txt",quote=F)