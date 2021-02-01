#Plot.figS3.AT.Snpeff.R
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



par(mfrow=c(1,4))


AT.all.mutations.annotation <- read.csv("../M2.mutations.full.list.annotation.csv",header=T)


sample.vec <- c(sort(unique(AT.all.mutations.annotation$Sample1)))
AT.all.mutations.annotation$Sample1 <- factor(AT.all.mutations.annotation$Sample1, levels=sample.vec)
AT.all.mutations.annotation$Sample2 <- factor(AT.all.mutations.annotation$Sample2, levels=sample.vec)
AT.all.mutations.annotation$Sample3 <- factor(AT.all.mutations.annotation$Sample3, levels=sample.vec)

AT.all.family <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Sample2!="NA")

AT.all.sbs <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Type=="SBS")
AT.all.indel <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Type!="SBS")
AT.all.insertion <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Type=="Insertion")
AT.all.deletion <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Type=="Deletion")

AT.all.HIGH <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Annotation_impact=="HIGH")
AT.all.MODERATE <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Annotation_impact=="MODERATE")
AT.all.LOW <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Annotation_impact=="LOW")
AT.all.MODIFIER <- subset(AT.all.mutations.annotation, AT.all.mutations.annotation$Annotation_impact=="MODIFIER")

No.mutations.per.sample1 <- tapply(AT.all.mutations.annotation$Chr, AT.all.mutations.annotation$Sample1, length)
No.mutations.per.sample2 <- tapply(AT.all.mutations.annotation$Chr, AT.all.mutations.annotation$Sample2, length)
No.mutations.per.sample3 <- tapply(AT.all.mutations.annotation$Chr, AT.all.mutations.annotation$Sample3, length)
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

No.homo.mutations.per.sample1 <- tapply(AT.all.mutations.annotation$Chr[AT.all.mutations.annotation$Zygosity1=="homo"], AT.all.mutations.annotation$Sample1[AT.all.mutations.annotation$Zygosity1=="homo"], length)
No.homo.mutations.per.sample2 <- tapply(AT.all.mutations.annotation$Chr[AT.all.mutations.annotation$Zygosity2=="homo"], AT.all.mutations.annotation$Sample2[AT.all.mutations.annotation$Zygosity2=="homo"], length)
No.homo.mutations.per.sample3 <- tapply(AT.all.mutations.annotation$Chr[AT.all.mutations.annotation$Zygosity3=="homo"], AT.all.mutations.annotation$Sample3[AT.all.mutations.annotation$Zygosity3=="homo"], length)
No.homo.mutations.per.sample1[is.na(No.homo.mutations.per.sample1)] <- 0
No.homo.mutations.per.sample2[is.na(No.homo.mutations.per.sample2)] <- 0
No.homo.mutations.per.sample3[is.na(No.homo.mutations.per.sample3)] <- 0
No.homo.mutations.per.sample <- No.homo.mutations.per.sample1 + No.homo.mutations.per.sample2 + No.homo.mutations.per.sample3

No.hetero.mutations.per.sample1 <- tapply(AT.all.mutations.annotation$Chr[AT.all.mutations.annotation$Zygosity1=="hetero"], AT.all.mutations.annotation$Sample1[AT.all.mutations.annotation$Zygosity1=="hetero"], length)
No.hetero.mutations.per.sample2 <- tapply(AT.all.mutations.annotation$Chr[AT.all.mutations.annotation$Zygosity2=="hetero"], AT.all.mutations.annotation$Sample2[AT.all.mutations.annotation$Zygosity2=="hetero"], length)
No.hetero.mutations.per.sample3 <- tapply(AT.all.mutations.annotation$Chr[AT.all.mutations.annotation$Zygosity3=="hetero"], AT.all.mutations.annotation$Sample3[AT.all.mutations.annotation$Zygosity3=="hetero"], length)
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

No.HIGH.per.sample1 <- tapply(AT.all.HIGH$Chr, AT.all.HIGH$Sample1, length)
No.HIGH.per.sample2 <- tapply(AT.all.HIGH$Chr, AT.all.HIGH$Sample2, length)
No.HIGH.per.sample3 <- tapply(AT.all.HIGH$Chr, AT.all.HIGH$Sample3, length)
No.HIGH.per.sample1[is.na(No.HIGH.per.sample1)] <- 0
No.HIGH.per.sample2[is.na(No.HIGH.per.sample2)] <- 0
No.HIGH.per.sample3[is.na(No.HIGH.per.sample3)] <- 0
No.HIGH.per.sample <- No.HIGH.per.sample1 + No.HIGH.per.sample2 + No.HIGH.per.sample3

No.MODERATE.per.sample1 <- tapply(AT.all.MODERATE$Chr, AT.all.MODERATE$Sample1, length)
No.MODERATE.per.sample2 <- tapply(AT.all.MODERATE$Chr, AT.all.MODERATE$Sample2, length)
No.MODERATE.per.sample3 <- tapply(AT.all.MODERATE$Chr, AT.all.MODERATE$Sample3, length)
No.MODERATE.per.sample1[is.na(No.MODERATE.per.sample1)] <- 0
No.MODERATE.per.sample2[is.na(No.MODERATE.per.sample2)] <- 0
No.MODERATE.per.sample3[is.na(No.MODERATE.per.sample3)] <- 0
No.MODERATE.per.sample <- No.MODERATE.per.sample1 + No.MODERATE.per.sample2 + No.MODERATE.per.sample3

No.LOW.per.sample1 <- tapply(AT.all.LOW$Chr, AT.all.LOW$Sample1, length)
No.LOW.per.sample2 <- tapply(AT.all.LOW$Chr, AT.all.LOW$Sample2, length)
No.LOW.per.sample3 <- tapply(AT.all.LOW$Chr, AT.all.LOW$Sample3, length)
No.LOW.per.sample1[is.na(No.LOW.per.sample1)] <- 0
No.LOW.per.sample2[is.na(No.LOW.per.sample2)] <- 0
No.LOW.per.sample3[is.na(No.LOW.per.sample3)] <- 0
No.LOW.per.sample <- No.LOW.per.sample1 + No.LOW.per.sample2 + No.LOW.per.sample3

No.MODIFIER.per.sample1 <- tapply(AT.all.MODIFIER$Chr, AT.all.MODIFIER$Sample1, length)
No.MODIFIER.per.sample2 <- tapply(AT.all.MODIFIER$Chr, AT.all.MODIFIER$Sample2, length)
No.MODIFIER.per.sample3 <- tapply(AT.all.MODIFIER$Chr, AT.all.MODIFIER$Sample3, length)
No.MODIFIER.per.sample1[is.na(No.MODIFIER.per.sample1)] <- 0
No.MODIFIER.per.sample2[is.na(No.MODIFIER.per.sample2)] <- 0
No.MODIFIER.per.sample3[is.na(No.MODIFIER.per.sample3)] <- 0
No.MODIFIER.per.sample <- No.MODIFIER.per.sample1 + No.MODIFIER.per.sample2 + No.MODIFIER.per.sample3



treat <- c(rep("Control",length=9),rep("Low",length=9),rep("Middle",length=9),rep("High",length=9))
gray <- c(rep(0,length=9),rep(0.4,length=9),rep(1.4,length=9),rep(2.0,length=9))
Group <- c(rep("Group1",18),rep("Group2",18))
Group <- factor(Group, levels=c("Group1", "Group2"))
Accumurate.gray <- gray*60
Accumurate.gray.revised <- c(rep(0,length=9),rep(23,length=9),rep(80,length=9),rep(114,length=9))
family <-c(rep("A01",length=3),rep("A02",length=3),rep("A03",length=3),rep("A11",length=3),rep("A12",length=3),rep("A13",length=3),rep("A21",length=3),rep("A22",length=3),rep("A23",length=3),rep("A31",length=3),rep("A32",length=3),rep("A33",length=3))
family <- factor(family, levels=c("A01","A02","A03","A11","A12","A13","A21","A22","A23","A31","A32","A33"))

mutation.annotation.count.frame <- data.frame(SampleID = sample.vec, Family = family,
	Treat = treat, Gray = gray, TotalGray = Accumurate.gray,
	Mutation.Count =  No.mutations.per.sample,
	Mutation.family.Count = No.fammily.mutations.per.sample,
	Mutation.homo.Count =  No.homo.mutations.per.sample,	
	Mutation.hetero.Count =  No.hetero.mutations.per.sample,	
	SBS.Count = No.sbs.per.sample,
	INDEL.Count = No.indel.per.sample,
	Insertion.Count = No.insertion.per.sample,
	Deletion.Count = No.deletion.per.sample,
	HIGH.Count = No.HIGH.per.sample,
	MODERATE.Count = No.MODERATE.per.sample,
	LOW.Count = No.LOW.per.sample,
	MODIFIER.Count = No.MODIFIER.per.sample,
	kakudai.group = Group
)

write.csv(mutation.annotation.count.frame, "No.mutations.annotation.summary.csv",quote=F, row.names=F)




#glm for HIGH impact
print("-------------------------------")
print("glm for HIGH")
out.HIGH.2nd.nls <- nls(mutation.annotation.count.frame$HIGH.Count ~　a+b*mutation.annotation.count.frame$Gray+c*(mutation.annotation.count.frame$Gray)^2, data= mutation.annotation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.HIGH.2nd.nls))
print("AIC for Second polymoninal model")
print(AIC(out.HIGH.2nd.nls))
out.glm.HIGH.gy <- glm(mutation.annotation.count.frame$HIGH.Count ~ mutation.annotation.count.frame$Gray, family = poisson)
print(summary(out.glm.HIGH.gy))
out.glm.HIGH.gy.quasi <- glm(mutation.annotation.count.frame$HIGH.Count ~ mutation.annotation.count.frame$Gray, family = quasipoisson)
print(summary(out.glm.HIGH.gy.quasi))
overdispertion.HIGH <- dispersiontest(out.glm.HIGH.gy)
print(overdispertion.HIGH)
out.glm.HIGH.gy.nb <- glm.nb(mutation.annotation.count.frame$HIGH.Count ~ mutation.annotation.count.frame$Gray)
print(summary(out.glm.HIGH.gy.nb))
odTest(out.glm.HIGH.gy.nb)
out.glmer.HIGH.gy <- glmer(HIGH.Count ~ Gray + (1|Family), family=poisson, data=mutation.annotation.count.frame)
print(summary(out.glmer.HIGH.gy))
out.glmer.HIGH.gy.nb <- glmer.nb(HIGH.Count ~ Gray + (1|Family), data=mutation.annotation.count.frame)
print(summary(out.glmer.HIGH.gy.nb))
print("")

print("-------------------------------")
print("glm for MODERATE")
out.MODERATE.2nd.nls <- nls(mutation.annotation.count.frame$MODERATE.Count ~　a+b*mutation.annotation.count.frame$Gray+c*(mutation.annotation.count.frame$Gray)^2, data= mutation.annotation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.MODERATE.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.MODERATE.2nd.nls))
out.glm.MODERATE.gy <- glm(mutation.annotation.count.frame$MODERATE.Count ~ mutation.annotation.count.frame$Gray, family = poisson)
print(summary(out.glm.MODERATE.gy))
out.glm.MODERATE.gy.quasi <- glm(mutation.annotation.count.frame$MODERATE.Count ~ mutation.annotation.count.frame$Gray, family = quasipoisson)
print(summary(out.glm.MODERATE.gy.quasi))
overdispertion.MODERATE <- dispersiontest(out.glm.MODERATE.gy)
print(overdispertion.MODERATE)
out.glm.MODERATE.gy.nb <- glm.nb(mutation.annotation.count.frame$MODERATE.Count ~ mutation.annotation.count.frame$Gray)
print(summary(out.glm.MODERATE.gy.nb))
odTest(out.glm.MODERATE.gy.nb)
out.glmer.MODERATE.gy <- glmer(MODERATE.Count ~ Gray + (1|Family), family=poisson, data=mutation.annotation.count.frame)
print(summary(out.glmer.MODERATE.gy))
out.glmer.MODERATE.gy.nb <- glmer.nb(MODERATE.Count ~ Gray + (1|Family), data=mutation.annotation.count.frame)
print(summary(out.glmer.MODERATE.gy.nb))
print("")

print("-------------------------------")
print("glm for LOW")
out.LOW.2nd.nls <- nls(mutation.annotation.count.frame$LOW.Count ~　a+b*mutation.annotation.count.frame$Gray+c*(mutation.annotation.count.frame$Gray)^2, data= mutation.annotation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.LOW.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.LOW.2nd.nls))
out.glm.LOW.gy <- glm(mutation.annotation.count.frame$LOW.Count ~ mutation.annotation.count.frame$Gray, family = poisson)
print(summary(out.glm.LOW.gy))
out.glm.LOW.gy.quasi <- glm(mutation.annotation.count.frame$LOW.Count ~ mutation.annotation.count.frame$Gray, family = quasipoisson)
print(summary(out.glm.LOW.gy.quasi))
overdispertion.LOW <- dispersiontest(out.glm.LOW.gy)
print(overdispertion.LOW)
out.glm.LOW.gy.nb <- glm.nb(mutation.annotation.count.frame$LOW.Count ~ mutation.annotation.count.frame$Gray)
print(summary(out.glm.LOW.gy.nb))
odTest(out.glm.LOW.gy.nb)
out.glmer.LOW.gy <- glmer(LOW.Count ~ Gray + (1|Family), family=poisson, data=mutation.annotation.count.frame)
print(summary(out.glmer.LOW.gy))
out.glmer.LOW.gy.nb <- glmer.nb(LOW.Count ~ Gray + (1|Family), data=mutation.annotation.count.frame)
print(summary(out.glmer.LOW.gy.nb))
print("")

print("-------------------------------")
print("glm for MODIFIER")
out.MODIFIER.2nd.nls <- nls(mutation.annotation.count.frame$MODIFIER.Count ~　a+b*mutation.annotation.count.frame$Gray+c*(mutation.annotation.count.frame$Gray)^2, data= mutation.annotation.count.frame, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.MODIFIER.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.MODIFIER.2nd.nls))
out.glm.MODIFIER.gy <- glm(mutation.annotation.count.frame$MODIFIER.Count ~ mutation.annotation.count.frame$Gray, family = poisson)
print(summary(out.glm.MODIFIER.gy))
out.glm.MODIFIER.gy.quasi <- glm(mutation.annotation.count.frame$MODIFIER.Count ~ mutation.annotation.count.frame$Gray, family = quasipoisson)
print(summary(out.glm.MODIFIER.gy.quasi))
overdispertion.MODIFIER <- dispersiontest(out.glm.MODIFIER.gy)
print(overdispertion.MODIFIER)
out.glm.MODIFIER.gy.nb <- glm.nb(mutation.annotation.count.frame$MODIFIER.Count ~ mutation.annotation.count.frame$Gray)
print(summary(out.glm.MODIFIER.gy.nb))
odTest(out.glm.MODIFIER.gy.nb)
out.glmer.MODIFIER.gy <- glmer(MODIFIER.Count ~ Gray + (1|Family), family=poisson, data=mutation.annotation.count.frame)
print(summary(out.glmer.MODIFIER.gy))
out.glmer.MODIFIER.gy.nb <- glmer.nb(MODIFIER.Count ~ Gray + (1|Family), data=mutation.annotation.count.frame)
print(summary(out.glmer.MODIFIER.gy.nb))
print("")



snp.indel.color <- "mediumpurple"
snp.color <- "deeppink3"
indel.color <- "mediumblue"
black.color <- "black"
samon.color <- rgb(1,0.75,0.85)

grid.newpage() #空の画面を作る
pushViewport(viewport(layout=grid.layout(2, 2))) #画面を区切る（今回は2行2列の4分割）


gg.plot.high <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
		axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
	geom_count(data=mutation.annotation.count.frame,aes(x=Gray, y= HIGH.Count), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.annotation.count.frame,aes(x=Gray, y= HIGH.Count), method = "glm", method.args = list(family = poisson), se =FALSE, col="black", lwt= 0.5) + ggtitle("a) High-impact mutations")

gg.plot.moderate <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
		axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
	geom_count(data=mutation.annotation.count.frame,aes(x=Gray, y= MODERATE.Count), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.annotation.count.frame,aes(x=Gray, y= MODERATE.Count), method = "glm", method.args = list(family = poisson), se =FALSE, col="black", lwt= 0.5) + ggtitle("b) Moderate-impact mutations")

gg.plot.low <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
		axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
	geom_count(data=mutation.annotation.count.frame,aes(x=Gray, y= LOW.Count), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.annotation.count.frame,aes(x=Gray, y= LOW.Count), method = "glm", method.args = list(family = poisson), se =FALSE, col="black", lwt= 0.5) + ggtitle("c) Low-impact mutations")

gg.plot.modifier <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
		axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
	geom_count(data=mutation.annotation.count.frame,aes(x=Gray, y= MODIFIER.Count), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=mutation.annotation.count.frame,aes(x=Gray, y= MODIFIER.Count), method = MASS::glm.nb, se =FALSE, col="black", lwt= 0.5) + ggtitle("d) Modifier-impact mutations")


print(gg.plot.high, vp=viewport(layout.pos.row=1, layout.pos.col=1)) #1行目の1列
print(gg.plot.moderate, vp=viewport(layout.pos.row=1, layout.pos.col=2)  ) #1行目の2列
print(gg.plot.low, vp=viewport(layout.pos.row=2, layout.pos.col=1)) #1行目の1列
print(gg.plot.modifier, vp=viewport(layout.pos.row=2, layout.pos.col=2)  ) #1行目の2列

