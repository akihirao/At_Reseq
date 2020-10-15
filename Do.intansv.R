#Do.intansv.R
#2020,10.2

library(intansv)

reference_folder <- "/zfs/Arabidopsis/Reference_v1.1"
breakdancer_folder <- "/zfs/Arabidopsis/work/At_Reseq/breakdancer_out"
pindel_folder <- "/zfs/Arabidopsis/work/At_Reseq/pindel_out"

args <- commandArgs()
target.sample <- args[6]

#target.sample <- "AT01"

print(target.sample)

target.sample.breakdancer.out <- paste(target.sample, "breakdancer", "out", sep=".")
breakdancer.file.path <- paste(breakdancer_folder, target.sample.breakdancer.out, sep ="/")

breakdancer <- readBreakDancer(breakdancer.file.path)
str(breakdancer)


pindle.file.path <- paste(pindel_folder, target.sample, sep ="/")
pindel <- readPindel(pindle.file.path)
str(pindel)


sv.breakdancer.pindel <- methodsMerge(breakdancer, pindel)
str(sv.breakdancer.pindel)

#browser()

output.del <- paste(target.sample,"intansv","del","out", sep=".")
output.dup <- paste(target.sample,"intansv","dup","out", sep=".")
output.inv <- paste(target.sample,"intansv","inv","out", sep=".")

write.table(sv.breakdancer.pindel$del, output.del, quote=F, row.names=F)
write.table(sv.breakdancer.pindel$dup, output.dup, quote=F, row.names=F)
write.table(sv.breakdancer.pindel$inv, output.inv, quote=F, row.names=F)