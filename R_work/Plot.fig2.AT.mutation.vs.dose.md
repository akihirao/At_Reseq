-   [Statistical modeling of the effect of radiation on the number of
    each type of
    mutation](#statistical-modeling-of-the-effect-of-radiation-on-the-number-of-each-type-of-mutation)
    -   [Loading packages](#loading-packages)
    -   [Loading data set](#loading-data-set)
    -   [Data preprocessing](#data-preprocessing)
    -   [Statistical analyses of the number of each type of
        mutation](#statistical-analyses-of-the-number-of-each-type-of-mutation)

------------------------------------------------------------------------

## Statistical modeling of the effect of radiation on the number of each type of mutation

### Loading packages

``` r
library(ggplot2)
library(stringr)
library(MASS)
library(glmmML)
library(lme4)
library(tidyverse)
library(nlstools)
library(ggforce)
library(grid)
library(AER)
library(pscl)
```

### Loading data set

``` r
AT_all_mutations <- read_csv("../M2.mutations.full.list.csv")
```

### Data preprocessing

``` r
sample.vec <- sort(unique(c(AT_all_mutations$Sample1,AT_all_mutations$Sample2,AT_all_mutations$Sample3)))

AT_all_mutations$Sample1 <- factor(AT_all_mutations$Sample1, levels=sample.vec)
AT_all_mutations$Sample2 <- factor(AT_all_mutations$Sample2, levels=sample.vec)
AT_all_mutations$Sample3 <- factor(AT_all_mutations$Sample3, levels=sample.vec)

AT_all_sbs <- AT_all_mutations %>% filter(Type=="SBS")
AT_all_indel <- AT_all_mutations %>% filter(Type!="SBS")
AT_all_insertion <- AT_all_mutations %>% filter(Type=="Insertion")
AT_all_deletion <- AT_all_mutations %>% filter(Type=="Deletion")

AT_all_family <- AT_all_mutations %>% filter(Sample2!="NA")

No_mutations_per_sample <- tapply(rep(AT_all_mutations$Chr,3), c(AT_all_mutations$Sample1,AT_all_mutations$Sample2,AT_all_mutations$Sample3), length)
No_mutations_per_sample[is.na(No_mutations_per_sample)] <- 0

No_family_mutations_per_sample <- tapply(rep(AT_all_family$Chr,3), c(AT_all_family$Sample1,AT_all_family$Sample2,AT_all_family$Sample3), length)
No_family_mutations_per_sample[is.na(No_family_mutations_per_sample)] <- 0

No_homo_mutations_per_sample <- tapply(c(AT_all_mutations$Chr[AT_all_mutations$Zygosity1=="homo"],AT_all_mutations$Chr[AT_all_mutations$Zygosity2=="homo"],AT_all_mutations$Chr[AT_all_mutations$Zygosity3=="homo"]), c(AT_all_mutations$Sample1[AT_all_mutations$Zygosity1=="homo"],AT_all_mutations$Sample2[AT_all_mutations$Zygosity2=="homo"],AT_all_mutations$Sample3[AT_all_mutations$Zygosity3=="homo"]), length)
No_homo_mutations_per_sample[is.na(No_homo_mutations_per_sample)] <- 0

No_hetero_mutations_per_sample <- tapply(c(AT_all_mutations$Chr[AT_all_mutations$Zygosity1=="hetero"],AT_all_mutations$Chr[AT_all_mutations$Zygosity2=="hetero"],AT_all_mutations$Chr[AT_all_mutations$Zygosity3=="hetero"]), c(AT_all_mutations$Sample1[AT_all_mutations$Zygosity1=="hetero"],AT_all_mutations$Sample2[AT_all_mutations$Zygosity2=="hetero"],AT_all_mutations$Sample3[AT_all_mutations$Zygosity3=="hetero"]), length)
No_hetero_mutations_per_sample[is.na(No_hetero_mutations_per_sample)] <- 0

No_sbs_per_sample <- tapply(rep(AT_all_sbs$Chr,3), c(AT_all_sbs$Sample1,AT_all_sbs$Sample2,AT_all_sbs$Sample3), length)
No_sbs_per_sample[is.na(No_sbs_per_sample)] <- 0

No_indel_per_sample <- tapply(rep(AT_all_indel$Chr,3), c(AT_all_indel$Sample1,AT_all_indel$Sample2,AT_all_indel$Sample3), length)
No_indel_per_sample[is.na(No_indel_per_sample)] <- 0

No_insertion_per_sample <- tapply(rep(AT_all_insertion$Chr,3), c(AT_all_insertion$Sample1,AT_all_insertion$Sample2,AT_all_insertion$Sample3), length)
No_insertion_per_sample[is.na(No_insertion_per_sample)] <- 0

No_deletion_per_sample <- tapply(rep(AT_all_deletion$Chr,3), c(AT_all_deletion$Sample1,AT_all_deletion$Sample2,AT_all_deletion$Sample3), length)
No_deletion_per_sample[is.na(No_deletion_per_sample)] <- 0

treat <- c(rep("Control",9),rep("Low",9),rep("Middle",9),rep("High",9))
gray <- c(rep(0,9),rep(0.4,9),rep(1.4,9),rep(2.0,9))
Accumurate.gray <- gray*60
Accumurate.gray.revised <- c(rep(0,9),rep(23,9),rep(80,9),rep(114,9))
family <-c(rep("A01",3),rep("A02",3),rep("A03",3),rep("A11",3),rep("A12",3),rep("A13",3),rep("A21",3),rep("A22",3),rep("A23",3),rep("A31",3),rep("A32",3),rep("A33",3))

mutation.count.frame <- data.frame(SampleID = sample.vec, Family = family,
    Treat = treat, Gray = gray, TotalGray = Accumurate.gray,
    Mutation.Count =  No_mutations_per_sample,
    Mutation.family.Count = No_family_mutations_per_sample,
    Mutation.homo.Count = No_homo_mutations_per_sample,
    Mutation.hetero.Count = No_hetero_mutations_per_sample,
    SBS.Count = No_sbs_per_sample,
    INDEL.Count = No_indel_per_sample,
    Insertion.Count = No_insertion_per_sample,
    Deletion.Count = No_deletion_per_sample
)
mutation.count.frame$Family <- factor(mutation.count.frame$Family, levels=c("A01","A02","A03","A11","A12","A13","A21","A22","A23","A31","A32","A33"))

write.csv(mutation.count.frame, "No.mutations.summary.csv",quote=F, row.names=F)
```

### Statistical analyses of the number of each type of mutation

#### Statistical modeling of the total number of mutations (SBSs + INDELs)

``` r
# GLM with a negative binominal distribution
glm.nb.total.mutation <- glm.nb(mutation.count.frame$Mutation.Count ~ mutation.count.frame$Gray)
print(summary(glm.nb.total.mutation))
```

    ## 
    ## Call:
    ## glm.nb(formula = mutation.count.frame$Mutation.Count ~ mutation.count.frame$Gray, 
    ##     init.theta = 10.95880274, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.0158  -0.9053  -0.1581   0.3501   3.2203  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                1.67664    0.12835   13.06   <2e-16 ***
    ## mutation.count.frame$Gray  0.95091    0.08989   10.58   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(10.9588) family taken to be 1)
    ## 
    ##     Null deviance: 162.742  on 35  degrees of freedom
    ## Residual deviance:  44.602  on 34  degrees of freedom
    ## AIC: 234.6
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  10.96 
    ##           Std. Err.:  5.56 
    ## 
    ##  2 x log-likelihood:  -228.602

``` r
# GLMM with a negative binominal distribution
glmer.nb.total.mutation <- glmer.nb(Mutation.Count ~ Gray + (1|Family), data=mutation.count.frame)
print(summary(glmer.nb.total.mutation))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: Negative Binomial(22.0047)  ( log )
    ## Formula: Mutation.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    230.7    237.0   -111.4    222.7       32 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.4383 -0.6915 -0.1216  0.3435  2.3673 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0.1071   0.3272  
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.5174     0.2057   7.376 1.63e-13 ***
    ## Gray          1.0472     0.1526   6.861 6.86e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.828

``` r
# GLM with a poisson distribution
glm.poisson.total.mutation <- glm(mutation.count.frame$Mutation.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(glm.poisson.total.mutation))
```

    ## 
    ## Call:
    ## glm(formula = mutation.count.frame$Mutation.Count ~ mutation.count.frame$Gray, 
    ##     family = poisson)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.5643  -1.4450  -0.2386   0.6267   4.5594  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                1.70805    0.09716   17.58   <2e-16 ***
    ## mutation.count.frame$Gray  0.92419    0.05986   15.44   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 374.180  on 35  degrees of freedom
    ## Residual deviance:  89.365  on 34  degrees of freedom
    ## AIC: 247.04
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
# GLMM with a poisson distribution
glmer.poisson.total.mutation <- glmer(Mutation.Count ~ Gray + (1|Family), family=poisson, data=mutation.count.frame)
print(summary(glmer.poisson.total.mutation))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: Mutation.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    236.2    241.0   -115.1    230.2       33 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2592 -0.8115 -0.1783  0.4220  3.7243 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0.1221   0.3495  
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.5135     0.2045   7.401 1.35e-13 ***
    ## Gray          1.0490     0.1520   6.899 5.23e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.826

``` r
# Non-lenear second polynominal model
nls.2nd.total.mutation <- nls(mutation.count.frame$Mutation.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
```

    ## 11176.99    (2.60e+00): par = (1 1 1)
    ## 1440.358    (3.98e-08): par = (4.803995 7.00905 3.728153)

``` r
print(summary(nls.2nd.total.mutation))
```

    ## 
    ## Formula: mutation.count.frame$Mutation.Count ~ a + b * mutation.count.frame$Gray + 
    ##     c * (mutation.count.frame$Gray)^2
    ## 
    ## Parameters:
    ##   Estimate Std. Error t value Pr(>|t|)  
    ## a    4.804      2.016   2.383   0.0231 *
    ## b    7.009      5.981   1.172   0.2496  
    ## c    3.728      2.924   1.275   0.2112  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 6.607 on 33 degrees of freedom
    ## 
    ## Number of iterations to convergence: 1 
    ## Achieved convergence tolerance: 3.984e-08

``` r
AIC(glm.nb.total.mutation, glmer.nb.total.mutation,glm.poisson.total.mutation, glmer.poisson.total.mutation, nls.2nd.total.mutation)
```

    ##                              df      AIC
    ## glm.nb.total.mutation         3 234.6023
    ## glmer.nb.total.mutation       4 230.7118
    ## glm.poisson.total.mutation    2 247.0440
    ## glmer.poisson.total.mutation  3 236.2454
    ## nls.2nd.total.mutation        4 242.9722

``` r
# Overdespersion test: poisson vs negative binomial
overdispertion.test.total.mutation <- dispersiontest(glm.poisson.total.mutation)
print(overdispertion.test.total.mutation)
```

    ## 
    ##  Overdispersion test
    ## 
    ## data:  glm.poisson.total.mutation
    ## z = 1.9002, p-value = 0.02871
    ## alternative hypothesis: true dispersion is greater than 1
    ## sample estimates:
    ## dispersion 
    ##   2.689317

``` r
odTest(glm.nb.total.mutation)
```

    ## Likelihood ratio test of H0: Poisson, as restricted NB model:
    ## n.b., the distribution of the test-statistic under H0 is non-standard
    ## e.g., see help(odTest) for details/references
    ## 
    ## Critical value of test statistic at the alpha= 0.05 level: 2.7055 
    ## Chi-Square Test Statistic =  14.4418 p-value = 7.228e-05

#### Statistical modeling of the number of SBS mutations

``` r
# GLM with a negative binomial distribution
glm.nb.sbs <- glm.nb(mutation.count.frame$SBS.Count ~ mutation.count.frame$Gray)
print(summary(glm.nb.sbs))
```

    ## 
    ## Call:
    ## glm.nb(formula = mutation.count.frame$SBS.Count ~ mutation.count.frame$Gray, 
    ##     init.theta = 16.38745177, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.5658  -0.8268  -0.1452   0.5906   2.5720  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                1.29349    0.13756   9.403   <2e-16 ***
    ## mutation.count.frame$Gray  0.95143    0.09181  10.363   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(16.3875) family taken to be 1)
    ## 
    ##     Null deviance: 169.058  on 35  degrees of freedom
    ## Residual deviance:  48.178  on 34  degrees of freedom
    ## AIC: 209.36
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  16.4 
    ##           Std. Err.:  11.3 
    ## 
    ##  2 x log-likelihood:  -203.364

``` r
# GLMM with a negative binomial distribution
glmer.nb.sbs <- glmer.nb(SBS.Count ~ Gray + (1|Family), data=mutation.count.frame)
print(summary(glmer.nb.sbs))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: Negative Binomial(23.7309)  ( log )
    ## Formula: SBS.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    209.3    215.7   -100.7    201.3       32 
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.55964 -0.73504  0.07293  0.47695  1.93402 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0.06978  0.2642  
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.1931     0.2006   5.949 2.69e-09 ***
    ## Gray          1.0103     0.1426   7.086 1.38e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.852

``` r
# GLM with a poisson distribution
glm.poisson.sbs <- glm(mutation.count.frame$SBS.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(glm.poisson.sbs))
```

    ## 
    ## Call:
    ## glm(formula = mutation.count.frame$SBS.Count ~ mutation.count.frame$Gray, 
    ##     family = poisson)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.7090  -1.0836  -0.1842   0.6782   3.1011  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                1.30004    0.11869   10.95   <2e-16 ***
    ## mutation.count.frame$Gray  0.94528    0.07284   12.98   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 274.278  on 35  degrees of freedom
    ## Residual deviance:  71.383  on 34  degrees of freedom
    ## AIC: 212.78
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
# GLMM with a poisson distribution
glmer.poisson.sbs <- glmer(SBS.Count ~ Gray + (1|Family), family=poisson, data=mutation.count.frame)
print(summary(glmer.poisson.sbs))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: SBS.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    210.9    215.6   -102.4    204.9       33 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7616 -0.8211  0.1204  0.5949  2.6418 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0.08325  0.2885  
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   1.1887     0.1988   5.978 2.25e-09 ***
    ## Gray          1.0120     0.1417   7.140 9.33e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.852

``` r
# Non-lenear second polynominal model
nls.2nd.sbs <- nls(mutation.count.frame$SBS.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
```

    ## 4426.669    (2.16e+00): par = (1 1 1)
    ## 780.2383    (7.92e-08): par = (3.489888 2.984207 3.563983)

``` r
print(summary(nls.2nd.sbs))
```

    ## 
    ## Formula: mutation.count.frame$SBS.Count ~ a + b * mutation.count.frame$Gray + 
    ##     c * (mutation.count.frame$Gray)^2
    ## 
    ## Parameters:
    ##   Estimate Std. Error t value Pr(>|t|)  
    ## a    3.490      1.484   2.352   0.0248 *
    ## b    2.984      4.402   0.678   0.5025  
    ## c    3.564      2.152   1.656   0.1072  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 4.862 on 33 degrees of freedom
    ## 
    ## Number of iterations to convergence: 1 
    ## Achieved convergence tolerance: 7.919e-08

``` r
AIC(glm.nb.sbs, glmer.nb.sbs, glm.poisson.sbs, glmer.poisson.sbs, nls.2nd.sbs)
```

    ##                   df      AIC
    ## glm.nb.sbs         3 209.3644
    ## glmer.nb.sbs       4 209.3500
    ## glm.poisson.sbs    2 212.7762
    ## glmer.poisson.sbs  3 210.8814
    ## nls.2nd.sbs        4 220.9025

``` r
# Overdespersion test: poisson vs negative binomial
overdispertion.test.sbs <- dispersiontest(glm.poisson.sbs)
print(overdispertion.test.sbs)
```

    ## 
    ##  Overdispersion test
    ## 
    ## data:  glm.poisson.sbs
    ## z = 2.0131, p-value = 0.02205
    ## alternative hypothesis: true dispersion is greater than 1
    ## sample estimates:
    ## dispersion 
    ##   1.910816

``` r
odTest(glm.nb.sbs)
```

    ## Likelihood ratio test of H0: Poisson, as restricted NB model:
    ## n.b., the distribution of the test-statistic under H0 is non-standard
    ## e.g., see help(odTest) for details/references
    ## 
    ## Critical value of test statistic at the alpha= 0.05 level: 2.7055 
    ## Chi-Square Test Statistic =  5.4118 p-value = 0.01

#### Statistical modeling of the number of INDEL mutations

``` r
# GLM with a negative binomial distribution
glm.nb.indel <- glm.nb(mutation.count.frame$INDEL.Count ~ mutation.count.frame$Gray)
print(summary(glm.nb.indel))
```

    ## 
    ## Call:
    ## glm.nb(formula = mutation.count.frame$INDEL.Count ~ mutation.count.frame$Gray, 
    ##     init.theta = 8.547471484, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -2.4840  -0.6683  -0.1123   0.4623   2.8570  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 0.5840     0.1957   2.984  0.00284 ** 
    ## mutation.count.frame$Gray   0.9053     0.1308   6.924  4.4e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(8.5475) family taken to be 1)
    ## 
    ##     Null deviance: 92.862  on 35  degrees of freedom
    ## Residual deviance: 40.387  on 34  degrees of freedom
    ## AIC: 171.85
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  8.55 
    ##           Std. Err.:  5.81 
    ## 
    ##  2 x log-likelihood:  -165.851

``` r
# GLMM with a negative binomial distribution
glmer.nb.indel <- glmer.nb(INDEL.Count ~ Gray + (1|Family), data=mutation.count.frame)
print(summary(glmer.nb.indel))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: Negative Binomial(8.5473)  ( log )
    ## Formula: INDEL.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    173.9    180.2    -82.9    165.9       32 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7918 -0.6020 -0.1107  0.4892  4.0552 
    ## 
    ## Random effects:
    ##  Groups Name        Variance  Std.Dev. 
    ##  Family (Intercept) 2.512e-12 1.585e-06
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.5840     0.1961   2.979   0.0029 ** 
    ## Gray          0.9053     0.1321   6.852 7.27e-12 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.870
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
# GLM with a poisson distribution
glm.poisson.indel <- glm(mutation.count.frame$INDEL.Count ~ mutation.count.frame$Gray, family = poisson)
print(summary(glm.poisson.indel))
```

    ## 
    ## Call:
    ## glm(formula = mutation.count.frame$INDEL.Count ~ mutation.count.frame$Gray, 
    ##     family = poisson)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.2794  -0.8021  -0.1348   0.6678   3.4598  
    ## 
    ## Coefficients:
    ##                           Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)                 0.6154     0.1691   3.639 0.000274 ***
    ## mutation.count.frame$Gray   0.8792     0.1051   8.363  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 143.321  on 35  degrees of freedom
    ## Residual deviance:  61.135  on 34  degrees of freedom
    ## AIC: 174.55
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
# GLMM with a poisson distribution
glmer.poisson.indel <- glmer(INDEL.Count ~ Gray + (1|Family), family=poisson, data=mutation.count.frame)
print(summary(glmer.poisson.indel))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: INDEL.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    175.9    180.6    -84.9    169.9       33 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.3507 -0.6240 -0.1717  0.5309  3.7704 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0.03302  0.1817  
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.5673     0.2071   2.739  0.00616 ** 
    ## Gray          0.9065     0.1356   6.685 2.31e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.889

``` r
# Non-lenear second polynominal model
nls.2nd.indel <- nls(mutation.count.frame$INDEL.Count ~　a+b*mutation.count.frame$Gray+c*(mutation.count.frame$Gray)^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
```

    ## 487.3088    (7.54e-01): par = (1 1 1)
    ## 310.6672    (1.54e-08): par = (1.314107 4.024844 0.1641699)

``` r
print(summary(nls.2nd.indel))
```

    ## 
    ## Formula: mutation.count.frame$INDEL.Count ~ a + b * mutation.count.frame$Gray + 
    ##     c * (mutation.count.frame$Gray)^2
    ## 
    ## Parameters:
    ##   Estimate Std. Error t value Pr(>|t|)
    ## a   1.3141     0.9364   1.403    0.170
    ## b   4.0248     2.7777   1.449    0.157
    ## c   0.1642     1.3578   0.121    0.904
    ## 
    ## Residual standard error: 3.068 on 33 degrees of freedom
    ## 
    ## Number of iterations to convergence: 1 
    ## Achieved convergence tolerance: 1.544e-08

``` r
AIC(glm.nb.indel, glmer.nb.indel, glm.poisson.indel, glmer.poisson.indel, nls.2nd.indel)
```

    ##                     df      AIC
    ## glm.nb.indel         3 171.8513
    ## glmer.nb.indel       4 173.8513
    ## glm.poisson.indel    2 174.5525
    ## glmer.poisson.indel  3 175.8874
    ## nls.2nd.indel        4 187.7509

``` r
overdispertion.test.indel <- dispersiontest(glm.poisson.indel)
print(overdispertion.test.indel)
```

    ## 
    ##  Overdispersion test
    ## 
    ## data:  glm.poisson.indel
    ## z = 1.2923, p-value = 0.09813
    ## alternative hypothesis: true dispersion is greater than 1
    ## sample estimates:
    ## dispersion 
    ##   1.692805

``` r
odTest(glm.nb.indel)
```

    ## Likelihood ratio test of H0: Poisson, as restricted NB model:
    ## n.b., the distribution of the test-statistic under H0 is non-standard
    ## e.g., see help(odTest) for details/references
    ## 
    ## Critical value of test statistic at the alpha= 0.05 level: 2.7055 
    ## Chi-Square Test Statistic =  4.7013 p-value = 0.01507

#### Statistical modeling of the number of deletion mutations

``` r
# GLMM with a negative binomial distribution
glm.nb.deletion <- glm.nb(Deletion.Count ~ Gray, data=mutation.count.frame)
print(summary(glm.nb.deletion))
```

    ## 
    ## Call:
    ## glm.nb(formula = Deletion.Count ~ Gray, data = mutation.count.frame, 
    ##     init.theta = 5.870337474, link = log)
    ## 
    ## Deviance Residuals: 
    ##      Min        1Q    Median        3Q       Max  
    ## -2.21048  -1.06826  -0.02644   0.41680   2.65342  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.4159     0.2171   1.916   0.0554 .  
    ## Gray          0.9670     0.1459   6.627 3.43e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(5.8703) family taken to be 1)
    ## 
    ##     Null deviance: 91.580  on 35  degrees of freedom
    ## Residual deviance: 43.852  on 34  degrees of freedom
    ## AIC: 170.85
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  5.87 
    ##           Std. Err.:  3.59 
    ## 
    ##  2 x log-likelihood:  -164.847

``` r
# GLMM with a negative binomial distribution
glmer.nb.deletion <- glmer.nb(Deletion.Count ~ Gray + (1|Family), data=mutation.count.frame)
print(summary(glmer.nb.deletion))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: Negative Binomial(5.939)  ( log )
    ## Formula: Deletion.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    172.8    179.2    -82.4    164.8       32 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.5654 -0.8993 -0.0324  0.4218  3.7948 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0.003849 0.06204 
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.4104     0.2601   1.578    0.115    
    ## Gray          0.9696     0.1653   5.865 4.49e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.890

``` r
# GLM with a poisson distribution
glm.poisson.deletion <- glm(Deletion.Count ~ Gray, family = poisson, data=mutation.count.frame)
print(summary(glm.poisson.deletion))
```

    ## 
    ## Call:
    ## glm(formula = Deletion.Count ~ Gray, family = poisson, data = mutation.count.frame)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -3.1365  -1.5087   0.0087   0.6303   3.3384  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.4624     0.1810   2.555   0.0106 *  
    ## Gray          0.9286     0.1114   8.335   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 153.466  on 35  degrees of freedom
    ## Residual deviance:  70.314  on 34  degrees of freedom
    ## AIC: 175.69
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
# GLMM with a poisson distribution
glmer.poisson.deletion <- glmer(Deletion.Count ~ Gray + (1|Family), family=poisson, data=mutation.count.frame)
print(summary(glmer.poisson.deletion))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: Deletion.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##    176.3    181.1    -85.2    170.3       33 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1271 -1.1085 -0.1113  0.5658  3.3878 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0.06154  0.2481  
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)   0.3713     0.2469   1.504    0.133    
    ## Gray          0.9812     0.1638   5.992 2.08e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.889

``` r
# Non-lenear second polynominal model
nls.2nd.deletion <- nls(Deletion.Count ~　a+b*Gray+c*Gray^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
```

    ## 428.1488    (6.12e-01): par = (1 1 1)
    ## 311.4543    (8.68e-08): par = (1.039201 3.909301 0.1432585)

``` r
print(summary(nls.2nd.deletion))
```

    ## 
    ## Formula: Deletion.Count ~ a + b * Gray + c * Gray^2
    ## 
    ## Parameters:
    ##   Estimate Std. Error t value Pr(>|t|)
    ## a   1.0392     0.9375   1.108    0.276
    ## b   3.9093     2.7812   1.406    0.169
    ## c   0.1433     1.3596   0.105    0.917
    ## 
    ## Residual standard error: 3.072 on 33 degrees of freedom
    ## 
    ## Number of iterations to convergence: 1 
    ## Achieved convergence tolerance: 8.683e-08

``` r
AIC(glm.nb.deletion, glmer.nb.deletion, glm.poisson.deletion, glmer.poisson.deletion, nls.2nd.deletion)
```

    ##                        df      AIC
    ## glm.nb.deletion         3 170.8470
    ## glmer.nb.deletion       4 172.8455
    ## glm.poisson.deletion    2 175.6850
    ## glmer.poisson.deletion  3 176.3260
    ## nls.2nd.deletion        4 187.8420

``` r
overdispertion.test.deletion <- dispersiontest(glm.poisson.deletion)
print(overdispertion.test.deletion)
```

    ## 
    ##  Overdispersion test
    ## 
    ## data:  glm.poisson.deletion
    ## z = 1.6776, p-value = 0.04672
    ## alternative hypothesis: true dispersion is greater than 1
    ## sample estimates:
    ## dispersion 
    ##   1.856946

``` r
odTest(glm.nb.deletion)
```

    ## Likelihood ratio test of H0: Poisson, as restricted NB model:
    ## n.b., the distribution of the test-statistic under H0 is non-standard
    ## e.g., see help(odTest) for details/references
    ## 
    ## Critical value of test statistic at the alpha= 0.05 level: 2.7055 
    ## Chi-Square Test Statistic =  6.8381 p-value = 0.004462

#### Statistical modeling of the number of insertion mutations

``` r
# GLM with a negative binomial distribution
glm.nb.insertion <- glm.nb(Insertion.Count ~ Gray,data= mutation.count.frame)
print(summary(glm.nb.insertion))
```

    ## 
    ## Call:
    ## glm.nb(formula = Insertion.Count ~ Gray, data = mutation.count.frame, 
    ##     init.theta = 10753.46085, link = log)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.0895  -0.8459  -0.7450   0.6666   1.4304  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.2818     0.4716  -2.718  0.00656 **
    ## Gray          0.3800     0.3325   1.143  0.25307   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for Negative Binomial(10753.46) family taken to be 1)
    ## 
    ##     Null deviance: 29.036  on 35  degrees of freedom
    ## Residual deviance: 27.695  on 34  degrees of freedom
    ## AIC: 62.311
    ## 
    ## Number of Fisher Scoring iterations: 1
    ## 
    ## 
    ##               Theta:  10753 
    ##           Std. Err.:  276901 
    ## Warning while fitting theta: iteration limit reached 
    ## 
    ##  2 x log-likelihood:  -56.311

``` r
# GLMM with a negative binomial distribution
glmer.nb.insertion <- glmer.nb(Insertion.Count ~ Gray + (1|Family), data=mutation.count.frame)
print(summary(glmer.nb.insertion))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: Negative Binomial(28441.79)  ( log )
    ## Formula: Insertion.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##     64.3     70.6    -28.2     56.3       32 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7704 -0.5982 -0.5268  0.7674  1.8256 
    ## 
    ## Random effects:
    ##  Groups Name        Variance  Std.Dev. 
    ##  Family (Intercept) 3.671e-10 1.916e-05
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.2818     0.4641  -2.762  0.00574 **
    ## Gray          0.3800     0.3276   1.160  0.24602   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.832
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
# GLM with a poisson distribution
glm.poisson.insertion <- glm(Insertion.Count ~ Gray, family = poisson, data=mutation.count.frame)
print(summary(glm.poisson.insertion))
```

    ## 
    ## Call:
    ## glm(formula = Insertion.Count ~ Gray, family = poisson, data = mutation.count.frame)
    ## 
    ## Deviance Residuals: 
    ##     Min       1Q   Median       3Q      Max  
    ## -1.0895  -0.8459  -0.7450   0.6667   1.4305  
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.2818     0.4715  -2.718  0.00656 **
    ## Gray          0.3800     0.3325   1.143  0.25305   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for poisson family taken to be 1)
    ## 
    ##     Null deviance: 29.037  on 35  degrees of freedom
    ## Residual deviance: 27.696  on 34  degrees of freedom
    ## AIC: 60.31
    ## 
    ## Number of Fisher Scoring iterations: 5

``` r
# GLMM with a poisson distribution
glmer.poisson.insertion <- glmer(Insertion.Count ~ Gray + (1|Family), family=poisson, data=mutation.count.frame)
print(summary(glmer.poisson.insertion))
```

    ## Generalized linear mixed model fit by maximum likelihood (Laplace
    ##   Approximation) [glmerMod]
    ##  Family: poisson  ( log )
    ## Formula: Insertion.Count ~ Gray + (1 | Family)
    ##    Data: mutation.count.frame
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##     62.3     67.1    -28.2     56.3       33 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -0.7704 -0.5982 -0.5268  0.7674  1.8256 
    ## 
    ## Random effects:
    ##  Groups Name        Variance Std.Dev.
    ##  Family (Intercept) 0        0       
    ## Number of obs: 36, groups:  Family, 12
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error z value Pr(>|z|)   
    ## (Intercept)  -1.2818     0.4716  -2.718  0.00656 **
    ## Gray          0.3800     0.3325   1.143  0.25306   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##      (Intr)
    ## Gray -0.837
    ## optimizer (Nelder_Mead) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
# Non-lenear second polynominal model
nls.2nd.insertion <- nls(Insertion.Count ~　a+b*Gray+c*Gray^2, data= mutation.count.frame, start=c(a=1, b=1, c=1),trace=T)
```

    ## 534.1488    (7.17e+00): par = (1 1 1)
    ## 10.18989    (5.07e-08): par = (0.2749064 0.1155432 0.0209113)

``` r
print(summary(nls.2nd.insertion))
```

    ## 
    ## Formula: Insertion.Count ~ a + b * Gray + c * Gray^2
    ## 
    ## Parameters:
    ##   Estimate Std. Error t value Pr(>|t|)
    ## a  0.27491    0.16958   1.621    0.115
    ## b  0.11554    0.50306   0.230    0.820
    ## c  0.02091    0.24592   0.085    0.933
    ## 
    ## Residual standard error: 0.5557 on 33 degrees of freedom
    ## 
    ## Number of iterations to convergence: 1 
    ## Achieved convergence tolerance: 5.068e-08

``` r
AIC(glm.nb.insertion, glmer.nb.insertion, glm.poisson.insertion, glmer.poisson.insertion, nls.2nd.insertion)
```

    ##                         df      AIC
    ## glm.nb.insertion         3 62.31060
    ## glmer.nb.insertion       4 64.31032
    ## glm.poisson.insertion    2 60.31015
    ## glmer.poisson.insertion  3 62.31015
    ## nls.2nd.insertion        4 64.72714

``` r
overdispertion.test.insertion <- dispersiontest(glm.poisson.insertion)
print(overdispertion.test.insertion)
```

    ## 
    ##  Overdispersion test
    ## 
    ## data:  glm.poisson.insertion
    ## z = -2.0204, p-value = 0.9783
    ## alternative hypothesis: true dispersion is greater than 1
    ## sample estimates:
    ## dispersion 
    ##  0.6769373

``` r
odTest(glm.nb.insertion)
```

    ## Likelihood ratio test of H0: Poisson, as restricted NB model:
    ## n.b., the distribution of the test-statistic under H0 is non-standard
    ## e.g., see help(odTest) for details/references
    ## 
    ## Critical value of test statistic at the alpha= 0.05 level: 2.7055 
    ## Chi-Square Test Statistic =  -4e-04 p-value = 0.5

#### Plotting the relationship between the number of each type of mutaion and radiation dose

``` r
sbs.indel.color <- "mediumpurple"
indel.color <- "mediumblue"
black.color <- "black"
samon.color <- rgb(1,0.75,0.85)

chairo.color <- rgb(128/255,64/255,0/255)
aka.color <- rgb(255/255,75/255,0/255)
total.color <- chairo.color

midori.color <- rgb(3/255,175/255,122/255)
pink.color <- rgb(255/255,128/255,130/255)
sbs.color <-pink.color
deletion.color <- ao.color <- rgb(0/255,90/255,255/255)
sora.color <- rgb(77/255,196/255,255/255)
insertion.color <- sora.color
orange.col <- rgb(246/255,170/255,0/255)
murasaki.col <- rgb(153/255,0/255,153/255)

gg.plot.total.mutation <- ggplot2::ggplot() +
    xlab("Radiation dose (Gy/d)") +
    ylab("No. mutations") +
    theme(axis.title.x = element_text(size=16), axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=16),axis.text.y = element_text(size=16)) +
    geom_count(data=mutation.count.frame,aes(x=Gray, y= Mutation.Count), pch=1, col=total.color) +
    geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Mutation.Count), method = lm, formula = y ~ exp(glmer.nb.total.mutation@beta[1] + glmer.nb.total.mutation@beta[2]*x), se =FALSE,col=total.color,lwd=0.8) + scale_size_continuous(breaks = seq(1,9,1)) +
    geom_count(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), pch=2, col=sbs.color) +　
    geom_smooth(data=mutation.count.frame,aes(x=Gray, y= SBS.Count), method = lm, formula = y ~ exp(glmer.nb.sbs@beta[1] + glmer.nb.sbs@beta[2]*x), se =FALSE,col=sbs.color, lwd=0.6) +
    geom_count(data=mutation.count.frame,aes(x=Gray, y= Insertion.Count), pch=3, col=insertion.color) +
    geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Insertion.Count), method = glm, method.args = list(family = "poisson"), se =FALSE, col=insertion.color,lwd=1,lty="dotted") +
    geom_count(data=mutation.count.frame,aes(x=Gray, y= Deletion.Count), pch=4, col=deletion.color) +
    geom_smooth(data=mutation.count.frame,aes(x=Gray, y= Deletion.Count), method = MASS::glm.nb, se =FALSE, col=deletion.color, lwd=0.6)

gg.plot.total.mutation
```

![](Plot.fig2.AT.mutation.vs.dose_files/figure-markdown_github/unnamed-chunk-39-1.png)
