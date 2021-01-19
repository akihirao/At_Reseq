#Plot.figS2.AT.Snpeff.R
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

AT.snpEff.count <- read.table("../snpEff_all_samples_summary_file.txt")



treat <- c(rep("Control",length=9),rep("Low",length=9),rep("Middle",length=9),rep("High",length=9))
gray <- c(rep(0,length=9),rep(0.4,length=9),rep(1.4,length=9),rep(2.0,length=9))
Group <- c(rep("Group1",18),rep("Group2",18))
Group <- factor(Group, levels=c("Group1", "Group2"))
Accumurate.gray <- gray*60


AT.snpEff.count.info <- data.frame(
	Sample=AT.snpEff.count$V1, 
	High.effect=AT.snpEff.count$V2, 
	Moderate.effect=AT.snpEff.count$V3, 
	Modifier.effect=AT.snpEff.count$V4,
	Low.effect=AT.snpEff.count$V5,
	Treat = treat, Gray = gray, TotalGray = Accumurate.gray
)



#glm for HIGH impact
print("-------------------------------")
print("glm for HIGH")
out.HIGH.2nd.nls <- nls(AT.snpEff.count.info$High.effect ~　a+b*AT.snpEff.count.info$Gray+c*(AT.snpEff.count.info$Gray)^2, data= AT.snpEff.count.info, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.HIGH.2nd.nls))
print("AIC for Second polymoninal model")
print(AIC(out.HIGH.2nd.nls))
out.glm.HIGH.gy <- glm(AT.snpEff.count.info$High.effect ~ AT.snpEff.count.info$Gray, family = poisson)
print(summary(out.glm.HIGH.gy))
out.glm.HIGH.gy.quasi <- glm(AT.snpEff.count.info$High.effect ~ AT.snpEff.count.info$Gray, family = quasipoisson)
print(summary(out.glm.HIGH.gy.quasi))
overdispertion.HIGH <- dispersiontest(out.glm.HIGH.gy)
print(overdispertion.HIGH)
out.glm.HIGH.gy.nb <- glm.nb(AT.snpEff.count.info$High.effect ~ AT.snpEff.count.info$Gray)
print(summary(out.glm.HIGH.gy.nb))
odTest(out.glm.HIGH.gy.nb)
out.glm.HIGH.zeroinfl.poisson <- zeroinfl(AT.snpEff.count.info$High.effect ~ AT.snpEff.count.info$Gray, dist="poisson")
print(summary(out.glm.HIGH.zeroinfl.poisson))
print("AIC for zefo-inflated Poisson model")
print(AIC(out.glm.HIGH.zeroinfl.poisson))
out.glm.HIGH.zeroinfl.nb <- zeroinfl(AT.snpEff.count.info$High.effect ~ AT.snpEff.count.info$Gray, dist="negbin")
print(summary(out.glm.HIGH.zeroinfl.nb))
print("AIC for zefo-inflated NB model")
print(AIC(out.glm.HIGH.zeroinfl.nb))
print("")

print("-------------------------------")
print("glm for MODERATE")
out.MODERATE.2nd.nls <- nls(AT.snpEff.count.info$Moderate.effect ~　a+b*AT.snpEff.count.info$Gray+c*(AT.snpEff.count.info$Gray)^2, data= AT.snpEff.count.info, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.MODERATE.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.MODERATE.2nd.nls))
out.glm.MODERATE.gy <- glm(AT.snpEff.count.info$Moderate.effect ~ AT.snpEff.count.info$Gray, family = poisson)
print(summary(out.glm.MODERATE.gy))
out.glm.MODERATE.gy.quasi <- glm(AT.snpEff.count.info$Moderate.effect ~ AT.snpEff.count.info$Gray, family = quasipoisson)
print(summary(out.glm.MODERATE.gy.quasi))
overdispertion.MODERATE <- dispersiontest(out.glm.MODERATE.gy)
print(overdispertion.MODERATE)
out.glm.MODERATE.gy.nb <- glm.nb(AT.snpEff.count.info$Moderate.effect ~ AT.snpEff.count.info$Gray)
print(summary(out.glm.MODERATE.gy.nb))
odTest(out.glm.MODERATE.gy.nb)
out.glm.MODERATE.zeroinfl.poisson <- zeroinfl(AT.snpEff.count.info$Moderate.effect ~ AT.snpEff.count.info$Gray, dist="poisson")
print(summary(out.glm.MODERATE.zeroinfl.poisson))
print("AIC for zefo-inflated Poisson model")
print(AIC(out.glm.MODERATE.zeroinfl.poisson))
out.glm.MODERATE.zeroinfl.nb <- zeroinfl(AT.snpEff.count.info$Moderate.effect ~ AT.snpEff.count.info$Gray, dist="negbin")
print(summary(out.glm.MODERATE.zeroinfl.nb))
print("AIC for zefo-inflated NB model")
print(AIC(out.glm.MODERATE.zeroinfl.nb))
print("")

print("-------------------------------")
print("glm for LOW")
out.LOW.2nd.nls <- nls(AT.snpEff.count.info$Low.effect ~　a+b*AT.snpEff.count.info$Gray+c*(AT.snpEff.count.info$Gray)^2, data= AT.snpEff.count.info, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.LOW.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.LOW.2nd.nls))
out.glm.LOW.gy <- glm(AT.snpEff.count.info$Low.effect ~ AT.snpEff.count.info$Gray, family = poisson)
print(summary(out.glm.LOW.gy))
out.glm.LOW.gy.quasi <- glm(AT.snpEff.count.info$Low.effect ~ AT.snpEff.count.info$Gray, family = quasipoisson)
print(summary(out.glm.LOW.gy.quasi))
overdispertion.LOW <- dispersiontest(out.glm.LOW.gy)
print(overdispertion.LOW)
out.glm.LOW.gy.nb <- glm.nb(AT.snpEff.count.info$Low.effect ~ AT.snpEff.count.info$Gray)
print(summary(out.glm.LOW.gy.nb))
odTest(out.glm.LOW.gy.nb)
out.glm.LOW.zeroinfl.poisson <- zeroinfl(AT.snpEff.count.info$Low.effect ~ AT.snpEff.count.info$Gray, dist="poisson")
print(summary(out.glm.LOW.zeroinfl.poisson))
print("AIC for zefo-inflated Poisson model")
print(AIC(out.glm.LOW.zeroinfl.poisson))
out.glm.LOW.zeroinfl.nb <- zeroinfl(AT.snpEff.count.info$Low.effect ~ AT.snpEff.count.info$Gray, dist="negbin")
print(summary(out.glm.LOW.zeroinfl.nb))
print("AIC for zefo-inflated NB model")
print(AIC(out.glm.LOW.zeroinfl.nb))
print("")

print("-------------------------------")
print("glm for MODIFIER")
out.MODIFIER.2nd.nls <- nls(AT.snpEff.count.info$Modifier.effect ~　a+b*AT.snpEff.count.info$Gray+c*(AT.snpEff.count.info$Gray)^2, data= AT.snpEff.count.info, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.MODIFIER.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.MODIFIER.2nd.nls))
out.glm.MODIFIER.gy <- glm(AT.snpEff.count.info$Modifier.effect ~ AT.snpEff.count.info$Gray, family = poisson)
print(summary(out.glm.MODIFIER.gy))
out.glm.MODIFIER.gy.quasi <- glm(AT.snpEff.count.info$Modifier.effect ~ AT.snpEff.count.info$Gray, family = quasipoisson)
print(summary(out.glm.MODIFIER.gy.quasi))
overdispertion.MODIFIER <- dispersiontest(out.glm.MODIFIER.gy)
print(overdispertion.MODIFIER)
out.glm.MODIFIER.gy.nb <- glm.nb(AT.snpEff.count.info$Modifier.effect ~ AT.snpEff.count.info$Gray)
print(summary(out.glm.MODIFIER.gy.nb))
odTest(out.glm.MODIFIER.gy.nb)
#out.glm.MODIFIER.zeroinfl.poisson <- zeroinfl(AT.snpEff.count.info$Modifier.effect ~ AT.snpEff.count.info$Gray, dist="poisson")
#print(summary(out.glm.MODIFIER.zeroinfl.poisson))
#print("AIC for zefo-inflated Poisson model")
#print(AIC(out.glm.MODIFIER.zeroinfl.poisson))
#out.glm.MODIFIER.zeroinfl.nb <- zeroinfl(AT.snpEff.count.info$Modifier.effect ~ AT.snpEff.count.info$Gray, dist="negbin")
#print(summary(out.glm.MODIFIER.zeroinfl.nb))
#print("AIC for zefo-inflated NB model")
#print(AIC(out.glm.MODIFIER.zeroinfl.nb))
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
	geom_count(data=AT.snpEff.count.info,aes(x=Gray, y= High.effect), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=AT.snpEff.count.info,aes(x=Gray, y= High.effect), method = "glm", method.args = list(family = poisson), se =TRUE, col="black", lwt= 0.5) + ggtitle("a) High impact")

gg.plot.moderate <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
		axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
	geom_count(data=AT.snpEff.count.info,aes(x=Gray, y= Moderate.effect), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=AT.snpEff.count.info,aes(x=Gray, y= Moderate.effect), method = MASS::glm.nb, se =TRUE, col="black", lwt= 0.5) + ggtitle("b) Moderate impact")

gg.plot.low <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
		axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
	geom_count(data=AT.snpEff.count.info,aes(x=Gray, y= Low.effect), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=AT.snpEff.count.info,aes(x=Gray, y= Low.effect), method = "glm", method.args = list(family = poisson), se =TRUE, col="black", lwt= 0.5) + ggtitle("c) Low impact")

gg.plot.modifier <- ggplot2::ggplot() + 
	xlab("Radiation dose (Gy/d)") + 
	ylab("No. mutations") +
	theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), 
		axis.text.x = element_text(size=14),axis.text.y = element_text(size=14)) +
	geom_count(data=AT.snpEff.count.info,aes(x=Gray, y= Modifier.effect), pch=16) +　scale_size_continuous(breaks = seq(1,9,1)) +
	geom_smooth(data=AT.snpEff.count.info,aes(x=Gray, y= Modifier.effect), method = MASS::glm.nb, se =TRUE, col="black", lwt= 0.5) + ggtitle("d) Modifier impact")


print(gg.plot.high, vp=viewport(layout.pos.row=1, layout.pos.col=1)) #1行目の1列
print(gg.plot.moderate, vp=viewport(layout.pos.row=1, layout.pos.col=2)  ) #1行目の2列
print(gg.plot.low, vp=viewport(layout.pos.row=2, layout.pos.col=1)) #1行目の1列
print(gg.plot.modifier, vp=viewport(layout.pos.row=2, layout.pos.col=2)  ) #1行目の2列

#plot(gg.plot.high)



