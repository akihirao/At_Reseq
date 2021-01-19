#Plot.fig2.AT.mutation.vs.dose.R
#by HIRAO Akira

library(ggplot2)
library(stringr)
library(MASS)
library(glmmML)
library(lme4)
library(tidyverse)
library("nlstools")
library(ggforce)
library(grid)
library(AER)
library(pscl) # testing overdispersion for glm.nb


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
No.fammily.mutations.per.sample <- No.family.mutations.per.sample1 + No.family.mutations.per.sample2 + No.family.mutations.per.sample3

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


treat <- c(rep("Control",length=9),rep("Low",length=9),rep("Middle",length=9),rep("High",length=9))
gray <- c(rep(0,length=9),rep(0.4,length=9),rep(0.8,length=9),rep(1.6,length=9))
gray.revised <- c(rep(0,length=9),rep(0.41,length=9),rep(1.4,length=9),rep(2.0,length=9))
Group <- c(rep("Group1",18),rep("Group2",18))
Group <- factor(Group, levels=c("Group1", "Group2"))
Accumurate.gray <- gray*60
Accumurate.gray.revised <- c(rep(0,length=9),rep(23,length=9),rep(80,length=9),rep(114,length=9))
family <-c(rep("A01",length=3),rep("A02",length=3),rep("A03",length=3),rep("A11",length=3),rep("A12",length=3),rep("A13",length=3),rep("A21",length=3),rep("A22",length=3),rep("A23",length=3),rep("A31",length=3),rep("A32",length=3),rep("A33",length=3))
family <- factor(family, levels=c("A01","A02","A03","A11","A12","A13","A21","A22","A23","A31","A32","A33"))

mutation.count.frame <- data.frame(SampleID = sample.vec, Family = family,
	Treat = treat, Gray = gray.revised, TotalGray = Accumurate.gray,
	Mutation.Count =  No.mutations.per.sample,
	Mutation.family.Count = No.fammily.mutations.per.sample,
	Mutation.homo.Count =  No.homo.mutations.per.sample,	
	Mutation.hetero.Count =  No.hetero.mutations.per.sample,	
	SBS.Count = No.sbs.per.sample,
	INDEL.Count = No.indel.per.sample,
	Insertion.Count = No.insertion.per.sample,
	Deletion.Count = No.deletion.per.sample,
	kakudai.group = Group
)


#glm for SBSs + INDELs
print("-------------------------------")
print("glm for total mutations")
out.total.mutation.2nd.nls <- nls(mutation.count.frame$Mutation.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.total.mutation.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.total.mutation.2nd.nls))
out.glm.total.mutation.gy <- glm(mutation.count.frame$Mutation.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(out.glm.total.mutation.gy))
out.glm.total.mutation.gy.quasi <- glm(mutation.count.frame$Mutation.Count ~ mutation.count.frame$Gray, family = quasipoisson)
print(summary(out.glm.total.mutation.gy.quasi))
overdispertion.test.total.mutation <- dispersiontest(out.glm.total.mutation.gy)
print(overdispertion.test.total.mutation)
out.glm.total.mutation.gy.nb <- glm.nb(mutation.count.frame$Mutation.Count ~ mutation.count.frame$Gray)
print(summary(out.glm.total.mutation.gy.nb))
odTest(out.glm.total.mutation.gy.nb)
out.glmer.total.mutation.gy.nb <- glmer.nb(Mutation.Count ~ Gray|Family, data=mutation.count.frame)
print(summary(out.glmer.total.mutation.gy.nb))
print("")

print("-------------------------------")
print("glm for SBS mutations")
out.sbs.2nd.nls <- nls(mutation.count.frame$SBS.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.sbs.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.sbs.2nd.nls))
out.glm.sbs.gy <- glm(mutation.count.frame$SBS.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(out.glm.sbs.gy))
overdispertion.test.sbs <- dispersiontest(out.glm.sbs.gy)
print(overdispertion.test.sbs)
out.glm.sbs.gy.nb <- glm.nb(mutation.count.frame$SBS.Count ~ mutation.count.frame$Gray)
print(summary(out.glm.sbs.gy.nb))
odTest(out.glm.sbs.gy.nb)
out.glmer.sbs.gy.nb <- glmer.nb(SBS.Count ~ Gray|Family, data=mutation.count.frame)
print(summary(out.glmer.sbs.gy.nb))
print("")

print("-------------------------------")
print("glm for INDEL mutations")
out.indel.2nd.nls <- nls(mutation.count.frame$INDEL.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.indel.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.indel.2nd.nls))
out.glm.indel.gy <- glm(mutation.count.frame$INDEL.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(out.glm.indel.gy))
overdispertion.test.indel <- dispersiontest(out.glm.indel.gy)
print(overdispertion.test.indel)
out.glm.indel.gy.nb <- glm.nb(mutation.count.frame$INDEL.Count ~ mutation.count.frame$Gray)
print(summary(out.glm.indel.gy.nb))
odTest(out.glm.indel.gy.nb)
out.glmer.indel.gy.nb <- glmer.nb(INDEL.Count ~ Gray|Family, data=mutation.count.frame)
print(summary(out.glmer.indel.gy.nb))
print("")

print("-------------------------------")
print("glm for insertion mutations")
out.insertion.2nd.nls <- nls(mutation.count.frame$Insertion.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.insertion.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.insertion.2nd.nls))
out.glm.insertion.gy <- glm(mutation.count.frame$Insertion.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(out.glm.insertion.gy))
overdispertion.test.insertion <- dispersiontest(out.glm.insertion.gy)
print(overdispertion.test.insertion)
out.glm.insertion.gy.nb <- glm.nb(mutation.count.frame$Insertion.Count ~ mutation.count.frame$Gray)
print(summary(out.glm.insertion.gy.nb))
odTest(out.glm.insertion.gy.nb)
out.glmer.insertion.gy.nb <- glmer.nb(Insertion.Count ~ Gray|Family, data=mutation.count.frame)
print(summary(out.glmer.insertion.gy.nb))
print("")

print("-------------------------------")
print("glm for deletion mutations")
out.deletion.2nd.nls <- nls(mutation.count.frame$Deletion.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.deletion.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.deletion.2nd.nls))
out.glm.deletion.gy <- glm(mutation.count.frame$Deletion.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(out.glm.deletion.gy))
overdispertion.test.deletion <- dispersiontest(out.glm.deletion.gy)
print(overdispertion.test.deletion)
out.glm.deletion.gy.nb <- glm.nb(mutation.count.frame$Deletion.Count ~ mutation.count.frame$Gray)
print(summary(out.glm.deletion.gy.nb))
odTest(out.glm.deletion.gy.nb)
out.glmer.deletion.gy.nb <- glmer.nb(Deletion.Count ~ Gray|Family, data=mutation.count.frame)
print(summary(out.glmer.deletion.gy.nb))
print("")


#grid.newpage() #空の画面を作る
#pushViewport(viewport(layout=grid.layout(1, 2))) #画面を区切る（今回は1行2列の2分割）

total.color <- "red"
sbs.indel.color <- "mediumpurple"
sbs.color <- "deeppink3"
indel.color <- "mediumblue"
black.color <- "black"
samon.color <- rgb(1,0.75,0.85)

gg.plot.total.mutation <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), 
		axis.text.x = element_text(size=16),axis.text.y = element_text(size=16)) +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= Mutation.Count), pch=1, col=total.color) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Mutation.Count), method = MASS::glm.nb, se =TRUE,col=total.color) + scale_size_continuous(breaks = seq(1,9,1)) +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), pch=2, col=sbs.color) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), method = MASS::glm.nb, se =TRUE, col=sbs.color,lwd=0.5) +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= Insertion.Count), pch=3, col="cyan4") +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Insertion.Count), method = MASS::glm.nb, se =FALSE, col="cyan4",lwd=0.5,lty="dotted") +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= Deletion.Count), pch=4, col="blue4") +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Deletion.Count), method = MASS::glm.nb, se =TRUE, col="blue4", lwd=0.5)
plot(gg.plot.total.mutation)


gg.plot.sbs.insertion.deletion <- ggplot2::ggplot() + 
	xlab("Radiation exposure (Gy/d)") + 
	ylab("No. mutations") +
#	ylim(0,15) + 
	theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16), 
		axis.text.x = element_text(size=16),axis.text.y = element_text(size=16)) +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), pch=2, col=sbs.color) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), method = MASS::glm.nb, se =TRUE, col=sbs.color,lwd=0.5) +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= Insertion.Count), pch=3, col="cyan4") +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Insertion.Count), method = MASS::glm.nb, se =FALSE, col="cyan4",lwd=0.5,lty="dotted") +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= Deletion.Count), pch=4, col="blue4") +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Deletion.Count), method = MASS::glm.nb, se =TRUE, col="blue4", lwd=0.5)
#plot(gg.plot.insertion.deletion)

#print(gg.plot.total.mutation, vp=viewport(layout.pos.row=1, layout.pos.col=1)) #1行目の1列
#print(gg.plot.sbs.insertion.deletion, vp=viewport(layout.pos.row=1, layout.pos.col=2)  ) #1行目の2列


gg.plot.snp.indel.kakudai <- ggplot2::ggplot() + 
	xlab("Radiation exposure (Gy/d)") + 
	ylab("No. mutations") +
	ylim(0,40) + 
	theme(axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), 
		axis.text.x = element_text(size=16),axis.text.y = element_text(size=16)) +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), pch=16, col=sbs.color) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), method = MASS::glm.nb, se =FALSE, col=sbs.color) +
	geom_count(data=mutation.count.frame,aes(x=Gray, y= INDEL.Count), pch=0, col=indel.color) +
	geom_smooth(data=mutation.count.frame,aes(x=Gray, y= INDEL.Count), method = MASS::glm.nb, se =FALSE, col=indel.color) + facet_zoom(xy = Group == "Group1") + xlim(0,0.4) + ylim(0,15) 
#plot(gg.plot.snp.indel.kakudai)


effective.genome.bp <- 117743122
spontanesous.mutarion.rate <- mutation.count.frame$Mutation.Count[1:9]/effective.genome.bp


