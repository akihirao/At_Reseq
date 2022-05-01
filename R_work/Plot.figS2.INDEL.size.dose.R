#AINDEL.size.dose.R
#2020,12.26

library(MASS)
library(glmmML)
library(lme4)
library(pscl)
library(AER)



#AT.all.mutations <- read.table("AT.all.mutations.txt",header=T)
AT.all.mutations <- read.csv("../M2.mutations.full.list.csv",header=T)
AT.all.mutations$Treatment <- factor(AT.all.mutations$Treatment, levels=c("Control","Low","Middle","High"))
sample.vec <- c(sort(unique(AT.all.mutations$Sample1)))
AT.all.mutations$Sample1 <- factor(AT.all.mutations$Sample1, levels=sample.vec)
AT.all.mutations$Sample2 <- factor(AT.all.mutations$Sample2, levels=sample.vec)
AT.all.mutations$Sample3 <- factor(AT.all.mutations$Sample3, levels=sample.vec)


AT.all.mutations$Type <- factor(AT.all.mutations$Type, levels=c("SBS","Deletion","Insertion"))
AT.all.mutations$Abs.Length <- abs(AT.all.mutations$Length)

AT.all.family <- subset(AT.all.mutations, AT.all.mutations$Sample2!="NA")

AT.all.sbs <- subset(AT.all.mutations, AT.all.mutations$Type=="SBS")
AT.all.indel <- subset(AT.all.mutations, AT.all.mutations$Type!="SBS")

AT.all.indel$Type <- factor(AT.all.indel$Type, levels=c("Deletion","Insertion"))
AT.all.indel$Abs.Length <- abs(AT.all.indel$Length)

AT.all.insertion <- subset(AT.all.indel, AT.all.indel$Type=="Insertion")
AT.all.deletion <- subset(AT.all.indel, AT.all.indel$Type=="Deletion")


sora.col <- rgb(77/255,196/255,255/255)
midori.col <- rgb(3/255,175/255,122/255)
orange.col <- rgb(246/255,170/255,0/255)
murasaki.col <- rgb(153/255,0/255,153/255)
akarui.murasaki.col <- rgb(201/255,172/255,230/255)
#col.parette <- c("#ff8082","#f6aa00","#03af7a","#005aff")
#col.parette <- c("#005aff","#03af7a","#f6aa00","red")
col.parette <- c(sora.col,midori.col,orange.col,murasaki.col)
#col.parette <- c("black","black","black","black")

#g.hist <- ggplot(AT.all.mutations, aes(x=Length, fill=Treatment)) #c("#ff8082","#f6aa00","#d8f255","#ffffff")
g.hist <- ggplot(AT.all.indel, aes(x=Length, fill=Treatment)) #c("#ff8082","#f6aa00","#d8f255","#ffffff")
#g2 <- g2 + geom_histogram(breaks=seq(-125,50,by=1), position = "dodge", alpha=0.8)#透明度80% alpha
#g.hist <- g.hist + geom_histogram(breaks=seq(-124.5,50.5,by=1), position = "identity",alpha=0.9) + xlab("INDEL size (bp)") + ylab("Frequency") +facet_grid(Treatment~.) + scale_fill_manual(values = col.parette) 
g.hist <- g.hist + geom_histogram(breaks=seq(-145.5,45.5,by=1), position = "identity",alpha=0.9) + xlab("Size of mutation (bp)") + ylab("Frequency") +facet_grid(Treatment~.) + scale_fill_manual(values = col.parette) 

plot(g.hist)




attach(AT.all.mutations)



out.size.2nd.nls <- nls(Abs.Length ~　a+b*Dose+c*Dose^2, data= AT.all.mutations, start=c(a=1, b=1, c=1),trace=T)
print(summary(out.size.2nd.nls))
print("AIC of Second polymoninal model")
print(AIC(out.size.2nd.nls))
out.glm.size.gy <- glm(Abs.Length~Dose, family = poisson)
print(summary(out.glm.size.gy))
overdispertion.test.size <- dispersiontest(out.glm.size.gy)
print(overdispertion.test.size)
out.glm.size.gy.nb <- glm.nb(Abs.Length~Dose)
print(summary(out.glm.size.gy.nb))
odTest(out.glm.size.gy.nb)
#out.glmer.size.gy.nb <- glmer.nb(SBS.Count ~ Gray|Family, data=mutation.count.frame)
#print(summary(out.glmer.sbs.gy.nb))
print("")



Abs.Length.hurdle.3.poisson <- hurdle(Abs.Length~Dose, dist="poisson")
print(summary(Abs.Length.hurdle.3.poisson))
print(AIC(Abs.Length.hurdle.3.poisson))

Abs.Length.hurdle.3.negbin <- hurdle(Abs.Length~Dose, dist="negbin")
print(summary(Abs.Length.hurdle.3.negbin))
print(AIC(Abs.Length.hurdle.3.poisson))


attach(AT.all.mutations)





