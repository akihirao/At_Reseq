#Plot.fig3.AT.homo2hetero.R
#by HIRAO Akira

library(stringr)
library(MASS)
library(glmmML)
library(tidyverse)
library(ggplot2)


library(RColorBrewer)

cols = brewer.pal(4, "Pastel1")  # preparing colar palette



AT.all.mutations <- read.csv("../M2.mutations.full.list.csv",header=T)


sample.vec <- c(sort(unique(AT.all.mutations$Sample1)))
AT.all.mutations$Sample1 <- factor(AT.all.mutations$Sample1, levels=sample.vec)
AT.all.mutations$Sample2 <- factor(AT.all.mutations$Sample2, levels=sample.vec)
AT.all.mutations$Sample3 <- factor(AT.all.mutations$Sample3, levels=sample.vec)

AT.all.family <- subset(AT.all.mutations, AT.all.mutations$Sample2!="NA")

AT.all.sbs <- subset(AT.all.mutations, AT.all.mutations$Type=="SBS")
AT.all.indel <- subset(AT.all.mutations, AT.all.mutations$Type!="SBS")
AT.all.insertion <- subset(AT.all.mutations, AT.all.mutations$Type=="Insertion")
AT.all.deletion <- subset(AT.all.mutations, AT.all.mutations$Type=="Deletion")

No.mutations.per.sample1 <- tapply(AT.all.mutations$Chr, AT.all.mutations$Sample1, length)
No.mutations.per.sample2 <- tapply(AT.all.mutations$Chr, AT.all.mutations$Sample2, length)
No.mutations.per.sample3 <- tapply(AT.all.mutations$Chr, AT.all.mutations$Sample3, length)
No.mutations.per.sample1[is.na(No.mutations.per.sample1)] <- 0
No.mutations.per.sample2[is.na(No.mutations.per.sample2)] <- 0
No.mutations.per.sample3[is.na(No.mutations.per.sample3)] <- 0
No.mutations.per.sample <- No.mutations.per.sample1 + No.mutations.per.sample2 + No.mutations.per.sample3

No.family.mutations.per.sample1 <- tapply(AT.all.family$Chr,AT.all.family$Sample1, length)
No.family.mutations.per.sample2 <- tapply(AT.all.family$Chr,AT.all.family$Sample2, length)
No.family.mutations.per.sample3 <- tapply(AT.all.family$Chr,AT.all.family$Sample3, length)
No.family.mutations.per.sample1[is.na(No.family.mutations.per.sample1)] <- 0
No.family.mutations.per.sample2[is.na(No.family.mutations.per.sample2)] <- 0
No.family.mutations.per.sample3[is.na(No.family.mutations.per.sample3)] <- 0
No.family.mutations.per.sample <- No.family.mutations.per.sample1 + No.family.mutations.per.sample2 + No.family.mutations.per.sample3

No.homo.mutations.per.sample1 <- tapply(AT.all.mutations$Chr[AT.all.mutations$Zygosity1=="homo"], AT.all.mutations$Sample1[AT.all.mutations$Zygosity1=="homo"], length)
No.homo.mutations.per.sample2 <- tapply(AT.all.mutations$Chr[AT.all.mutations$Zygosity2=="homo"], AT.all.mutations$Sample2[AT.all.mutations$Zygosity2=="homo"], length)
No.homo.mutations.per.sample3 <- tapply(AT.all.mutations$Chr[AT.all.mutations$Zygosity3=="homo"], AT.all.mutations$Sample3[AT.all.mutations$Zygosity3=="homo"], length)
No.homo.mutations.per.sample1[is.na(No.homo.mutations.per.sample1)] <- 0
No.homo.mutations.per.sample2[is.na(No.homo.mutations.per.sample2)] <- 0
No.homo.mutations.per.sample3[is.na(No.homo.mutations.per.sample3)] <- 0
No.homo.mutations.per.sample <- No.homo.mutations.per.sample1 + No.homo.mutations.per.sample2 + No.homo.mutations.per.sample3

No.hetero.mutations.per.sample1 <- tapply(AT.all.mutations$Chr[AT.all.mutations$Zygosity1=="hetero"], AT.all.mutations$Sample1[AT.all.mutations$Zygosity1=="hetero"], length)
No.hetero.mutations.per.sample2 <- tapply(AT.all.mutations$Chr[AT.all.mutations$Zygosity2=="hetero"], AT.all.mutations$Sample2[AT.all.mutations$Zygosity2=="hetero"], length)
No.hetero.mutations.per.sample3 <- tapply(AT.all.mutations$Chr[AT.all.mutations$Zygosity3=="hetero"], AT.all.mutations$Sample3[AT.all.mutations$Zygosity3=="hetero"], length)
No.hetero.mutations.per.sample1[is.na(No.hetero.mutations.per.sample1)] <- 0
No.hetero.mutations.per.sample2[is.na(No.hetero.mutations.per.sample2)] <- 0
No.hetero.mutations.per.sample3[is.na(No.hetero.mutations.per.sample3)] <- 0
No.hetero.mutations.per.sample <- No.hetero.mutations.per.sample1 + No.hetero.mutations.per.sample2 + No.hetero.mutations.per.sample3

No.sbs.per.sample1 <- tapply(AT.all.sbs$Chr, AT.all.sbs$Sample1, length)
No.sbs.per.sample2 <- tapply(AT.all.sbs$Chr, AT.all.sbs$Sample2, length)
No.sbs.per.sample3 <- tapply(AT.all.sbs$Chr, AT.all.sbs$Sample3, length)
No.sbs.per.sample1[is.na(No.sbs.per.sample1)] <- 0
No.sbs.per.sample2[is.na(No.sbs.per.sample2)] <- 0
No.sbs.per.sample3[is.na(No.sbs.per.sample3)] <- 0
No.sbs.per.sample <- No.sbs.per.sample1 + No.sbs.per.sample2 + No.sbs.per.sample3

No.indel.per.sample1 <- tapply(AT.all.indel$Chr, AT.all.indel$Sample1, length)
No.indel.per.sample2 <- tapply(AT.all.indel$Chr, AT.all.indel$Sample2, length)
No.indel.per.sample3 <- tapply(AT.all.indel$Chr, AT.all.indel$Sample3, length)
No.indel.per.sample1[is.na(No.indel.per.sample1)] <- 0
No.indel.per.sample2[is.na(No.indel.per.sample2)] <- 0
No.indel.per.sample3[is.na(No.indel.per.sample3)] <- 0
No.indel.per.sample <- No.indel.per.sample1 + No.indel.per.sample2 + No.indel.per.sample3

No.insertion.per.sample1 <- tapply(AT.all.insertion$Chr, AT.all.insertion$Sample1, length)
No.insertion.per.sample2 <- tapply(AT.all.insertion$Chr, AT.all.insertion$Sample2, length)
No.insertion.per.sample3 <- tapply(AT.all.insertion$Chr, AT.all.insertion$Sample3, length)
No.insertion.per.sample1[is.na(No.insertion.per.sample1)] <- 0
No.insertion.per.sample2[is.na(No.insertion.per.sample2)] <- 0
No.insertion.per.sample3[is.na(No.insertion.per.sample3)] <- 0
No.insertion.per.sample <- No.insertion.per.sample1 + No.insertion.per.sample2 + No.insertion.per.sample3

No.deletion.per.sample1 <- tapply(AT.all.deletion$Chr, AT.all.deletion$Sample1, length)
No.deletion.per.sample2 <- tapply(AT.all.deletion$Chr, AT.all.deletion$Sample2, length)
No.deletion.per.sample3 <- tapply(AT.all.deletion$Chr, AT.all.deletion$Sample3, length)
No.deletion.per.sample1[is.na(No.deletion.per.sample1)] <- 0
No.deletion.per.sample2[is.na(No.deletion.per.sample2)] <- 0
No.deletion.per.sample3[is.na(No.deletion.per.sample3)] <- 0
No.deletion.per.sample <- No.deletion.per.sample1 + No.deletion.per.sample2 + No.deletion.per.sample3

Treat.vec <- c("Control","Low","Middle","High")
treat <- c(rep("Control",length=9),rep("Low",length=9),rep("Middle",length=9),rep("High",length=9))
gray <- c(rep(0,length=9),rep(0.4,length=9),rep(1.4,length=9),rep(2.0,length=9))
Group <- c(rep("Group1",18),rep("Group2",18))
Group <- factor(Group, levels=c("Group1", "Group2"))
Accumurate.gray <- gray*60
Accumurate.gray.revised <- c(rep(0,length=9),rep(23,length=9),rep(80,length=9),rep(114,length=9))
family <-c(rep("A01",length=3),rep("A02",length=3),rep("A03",length=3),rep("A11",length=3),rep("A12",length=3),rep("A13",length=3),rep("A21",length=3),rep("A22",length=3),rep("A23",length=3),rep("A31",length=3),rep("A32",length=3),rep("A33",length=3))
family <- factor(family, levels=c("A01","A02","A03","A11","A12","A13","A21","A22","A23","A31","A32","A33"))

mutation.count.frame <- data.frame(SampleID = sample.vec, Family = family,
	Treat = treat, Gray = gray, TotalGray = Accumurate.gray,
	Mutation.Count =  No.mutations.per.sample,
	Mutation.family.Count = No.family.mutations.per.sample,
	Mutation.homo.Count =  No.homo.mutations.per.sample,	
	Mutation.hetero.Count =  No.hetero.mutations.per.sample,
	HomoHeteroRatio=No.homo.mutations.per.sample/No.hetero.mutations.per.sample,
	SBS.Count = No.sbs.per.sample,
	INDEL.Count = No.indel.per.sample,
	Insertion.Count = No.insertion.per.sample,
	Deletion.Count = No.deletion.per.sample,
	kakudai.group = Group
)
mutation.count.frame$Treat <- factor(mutation.count.frame$Treat, levels=Treat.vec)

#===========================


HomoHeteroRatio.treat.mean <- tapply(mutation.count.frame$HomoHeteroRatio, mutation.count.frame$Treat,mean)
mutation.total.mean <- tapply(mutation.count.frame$Mutation.Count, mutation.count.frame$Treat,mean)
mutation.total.sd <- tapply(mutation.count.frame$Mutation.Count, mutation.count.frame$Treat,sd)
mutation.homo.sum <- tapply(mutation.count.frame$Mutation.homo.Count, mutation.count.frame$Treat,sum)
mutation.homo.mean <- tapply(mutation.count.frame$Mutation.homo.Count, mutation.count.frame$Treat,mean)
mutation.homo.sd <- tapply(mutation.count.frame$Mutation.homo.Count, mutation.count.frame$Treat,sd)
mutation.hetero.sum <- tapply(mutation.count.frame$Mutation.hetero.Count, mutation.count.frame$Treat,sum)
mutation.hetero.mean <- tapply(mutation.count.frame$Mutation.hetero.Count, mutation.count.frame$Treat,mean)
mutation.hetero.sd <- tapply(mutation.count.frame$Mutation.hetero.Count, mutation.count.frame$Treat,sd)

kai.square.test.homo.hetero.control <- chisq.test(x=c(mutation.homo.sum[1],mutation.hetero.sum[1],p=c(1/3,2/3)))
print(kai.square.test.homo.hetero.control)
kai.square.test.homo.hetero.low <- chisq.test(x=c(mutation.homo.sum[2],mutation.hetero.sum[2],p=c(1/3,2/3)))
print(kai.square.test.homo.hetero.low)
kai.square.test.homo.hetero.middle <- chisq.test(x=c(mutation.homo.sum[3],mutation.hetero.sum[3],p=c(1/3,2/3)))
print(kai.square.test.homo.hetero.middle)
kai.square.test.homo.hetero.high <- chisq.test(x=c(mutation.homo.sum[4],mutation.hetero.sum[4],p=c(1/3,2/3)))
print(kai.square.test.homo.hetero.high)




Mutation.homo.hetero.total.mean <- c(mutation.homo.mean,mutation.hetero.mean, mutation.total.mean)
Mutation.homo.hetero.total.sd <- c(mutation.homo.sd,mutation.hetero.sd, mutation.total.sd)

Mutation.homo.hetero.mean <- c(mutation.homo.mean,mutation.hetero.mean)
Mutation.homo.hetero.sd <- c(mutation.homo.sd,mutation.hetero.sd)


Treat.lab <- c("Control: 0.0 Gy/day","Low: 0.4 Gy/day","Middle: 1.4 Gy/day","High: 2.0 Gy/day")
mutation.treat.mean.info <- data.frame(Mutation.mean= Mutation.homo.hetero.mean, Mutation.sd=Mutation.homo.hetero.sd, Treatment=rep(Treat.vec,2), Zygosity=c(rep("Homozygous",4),rep("Heterozygous",4)))
mutation.treat.mean.info$Treatment <- factor(mutation.treat.mean.info$Treatment, levels=Treat.vec)
mutation.treat.mean.info$Zygosity <- factor(mutation.treat.mean.info$Zygosity, levels=c("Homozygous","Heterozygous"))

#col.parette <- c("#ff8082","#f6aa00","#03af7a","#005aff")
col.parette <- c("#005aff","#03af7a","#f6aa00","red")
g <- ggplot(mutation.treat.mean.info, aes(x = Zygosity, y = Mutation.mean, fill = Treatment)) +
	theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), 
		axis.text.x = element_text(size=16),axis.text.y = element_text(size=16))	
g <- g + geom_bar(stat = "identity", position ="dodge") + geom_errorbar(aes(ymin = Mutation.mean - Mutation.sd, ymax = Mutation.mean + Mutation.sd, width=0.3), position = position_dodge(width = 0.9)) + ylab("Frequency") + scale_fill_manual(values = col.parette)

plot(g)



kruskal.test.homo.out <- kruskal.test(mutation.count.frame$Mutation.homo.Count~mutation.count.frame$Treat)
print(kruskal.test.homo.out)
pairwise.wilcox.test.homo <- pairwise.wilcox.test(mutation.count.frame$Mutation.homo.Count,mutation.count.frame$Treat,exact=F)
print(pairwise.wilcox.test.homo)

kruskal.test.hetero.out <- kruskal.test(mutation.count.frame$Mutation.hetero.Coun~mutation.count.frame$Treat)
print(kruskal.test.hetero.out)
pairwise.wilcox.test.hetero <- pairwise.wilcox.test(mutation.count.frame$Mutation.hetero.Coun,mutation.count.frame$Treat,exact=F)
print(pairwise.wilcox.test.hetero)



No.sample <- length(mutation.count.frame$Sample)
binom.test.p <- numeric(No.sample)
binom.test.CI.lower <- numeric(No.sample)
binom.test.CI.upper <- numeric(No.sample)

for(i in 1:No.sample){
	target.hetero.count <- mutation.count.frame$Mutation.hetero.Count[i]
	target.homo.count <- mutation.count.frame$Mutation.homo.Count[i]
	binom.test.out <- binom.test(target.hetero.count, target.hetero.count+target.homo.count, 2/3, alternative="greater")
	binom.test.p[i] <- binom.test.out[3]
	binom.test.CI.lower[i] <- binom.test.out[[4]][1]
	binom.test.CI.upper[i] <- binom.test.out[[4]][2]
}
binom.test.p.vec <- unlist(binom.test.p)

hetero.homo.ratio.out.frame <- data.frame(
	SampleID = sample.vec, Family = family,
	Treat = treat, Gray = gray,
	Mutation.Count =  No.mutations.per.sample,
	Mutation.homo.Count =  No.homo.mutations.per.sample,	
	Mutation.hetero.Count =  No.hetero.mutations.per.sample,
	Proportion.hetero=No.hetero.mutations.per.sample/No.mutations.per.sample,
	Binom.test.P=binom.test.p.vec,
	Binom.test.CI.lower=binom.test.CI.lower,
	Binom.test.CI.upper=binom.test.CI.upper
)

write.csv(hetero.homo.ratio.out.frame,"hetero.homo.test.summary.csv",quote=F, row.names=F)
