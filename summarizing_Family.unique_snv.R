#summarizing_Family.unique_snv.R
#2020,3.13 (the last day of Suppin)

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
sample13 <- args[16]
sample14 <- args[17]
sample15 <- args[18]
sample16 <- args[19]
sample17 <- args[20]
sample18 <- args[21]
sample19 <- args[22]
sample20 <- args[23]
sample21 <- args[24]
sample22 <- args[25]
sample23 <- args[26]
sample24 <- args[27]
sample25 <- args[28]
sample26 <- args[29]
sample27 <- args[30]
sample28 <- args[31]
sample29 <- args[32]
sample30 <- args[33]
sample31 <- args[34]
sample32 <- args[35]
sample33 <- args[36]
sample34 <- args[37]
sample35 <- args[38]
sample36 <- args[39]
sample37 <- args[40]
sample38 <- args[41]
sample39 <- args[42]
sample40 <- args[43]
sample41 <- args[44]
sample42 <- args[45]
sample43 <- args[46]
sample44 <- args[47]
sample45 <- args[48]
sample46 <- args[49]
sample47 <- args[50]
sample48 <- args[51]

input_snvs_file <- paste(category,".family.unique_snvs_list.txt",sep="")
family.snvs.list <- read.table(input_snvs_file)
no.total.snv <- length(family.snvs.list[,1])

family.snvs.list$V4 <- factor(family.snvs.list$V4, 
	levels=c(category, sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12, sample13, sample14, sample15, sample16, sample17, sample18, sample19, sample20, sample21, sample22, sample23, sample24, sample25, sample26, sample27, sample28, sample29, sample30, sample31, sample32, sample33, sample34, sample35, sample36, sample37, sample38, sample39, sample40, sample41, sample42, sample43, sample44, sample45, sample46, sample47, sample48))
family.no_unique_snvs <- table(family.snvs.list$V4)

#print(family.no_unique_snps)


input_snvs_file.filtered <- paste(category,".family.unique_snvs_list.filtered.txt",sep="")
family.snvs.list.filtered <- read.table(input_snvs_file.filtered)
no.total.snv.filtered <- length(family.snvs.list.filtered[,1])

family.snvs.list.filtered$V4 <- factor(family.snvs.list.filtered$V4, 
	levels=c(category, sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12, sample13, sample14, sample15, sample16, sample17, sample18, sample19, sample20, sample21, sample22, sample23, sample24, sample25, sample26, sample27, sample28, sample29, sample30, sample31, sample32, sample33, sample34, sample35, sample36, sample37, sample38, sample39, sample40, sample41, sample42, sample43, sample44, sample45, sample46, sample47, sample48))
family.no_unique_snvs.filtered <- table(family.snvs.list.filtered$V4)

#print(family.no_unique_snps.filtered)


#output.no_unique.snv <- rbind(family.no_unique_snps, family.no_unique_indels,family.no_unique_snps.filtered,family.no_unique_indels.filtered)
#row.names(output.no_unique.snv) <- c("no.unique.snps","no.unique.indels","no.unique.snps.filtered","no.unique.indels.filtered")

output.no_unique.snv <- rbind(family.no_unique_snvs, family.no_unique_snvs.filtered)
row.names(output.no_unique.snv) <- c("no.unique.snvs","no.unique.snvs.filtered")


print(output.no_unique.snv)
write.table(output.no_unique.snv, "summary.no_snv.txt",quote=F)