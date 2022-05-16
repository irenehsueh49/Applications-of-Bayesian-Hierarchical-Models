---
title: "Irene Hsueh's BS 849 Homework 5"
author: "Irene Hsueh"
date: "2/28/2022"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
library(rjags)
library(coda)
library(formatR)
set.seed(1234)
```


### Long Life Family Study Dataset 
```{r}
long_life <- read.csv("C:/Irene Hsueh's Documents/MS Applied Biostatistics/BS 849 - Bayesian Modeling for Biomedical Research & Public Health/Class 5 - Applications of Hierarchical Models/Lecture/llfs_data.csv") %>% 
#Excluding Observations with NA
  na.omit() %>% 
#Selecting Variables
  dplyr::select(id=pedid, 
                sex=sex, 
                pc1=pc1, 
                pc2=pc2, 
                outcome=outcome, 
                SNP1=rs1046914, 
                SNP2=rs187337462, 
                snp4=rs904957,
                SNP4=rs1129227, 
                SNP5=rs1129226, 
                SNP6=rs7091963,
                SNP7=rs3750681) %>%
#Recoding Sex Variable
  mutate(sex = recode(sex, "1"="0", "2"="1"))

head(long_life, 10)

#Family Indicator 
family_indicator <- rep(NA, length(long_life$id)) 
for(i in 1:length(unique(long_life$id)))
  {
  family_indicator[which(long_life$id == unique(long_life$id)[i])] <-
    rep(i, length(which(long_life$id == unique(long_life$id)[i]
                        )
                  )
        )
  }
table(table(family_indicator))


#Creating SNP Data
SNP_data <- list(N = nrow(long_life), 
                 n_family=length(unique(family_indicator)),
                 family=family_indicator,
                 sex=long_life$sex, 
                 pc1=long_life$pc1, 
                 pc2=long_life$pc2, 
                 outcome=long_life$outcome, 
                 snp1 = long_life$SNP1, 
                 snp2 = long_life$SNP2, 
                 snp3 = long_life$snp4, 
                 snp4 = long_life$SNP4, 
                 snp5 = long_life$SNP5, 
                 snp6 = long_life$SNP6, 
                 snp7 = long_life$SNP7)
```



### SNP1 Analysis
```{r}
snp1_model_bugs <- 
"model{
#Logistic Model
for (i in 1:N){
  outcome[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_sex*sex[i] + beta_pc1*(pc1[i]-pc1_mean) +  beta_pc2*(pc2[i]-pc2_mean) +  beta_snp1*(snp1[i]) + beta_family[family[i]]
}

#Random Effect Per Family
for (i in 1:n_family){
  beta_family[i] ~ dnorm(beta0, tau)
}

#Prior Distributions
beta0 ~ dnorm(0, 0.01)
beta_sex ~ dnorm(0, 0.01)
beta_pc1 ~ dnorm(0, 0.01)
beta_pc2 ~ dnorm(0, 0.01)
beta_snp1 ~ dnorm(0, 0.01)
tau ~ dgamma(1, 1)

OR1 <- exp(beta_snp1)
pc1_mean <- mean(pc1[])
pc2_mean <- mean(pc2[])
}"

snp1_model <- jags.model(textConnection(snp1_model_bugs), data=SNP_data, n.adapt=1000, n.chains=3)
snp1_model_gibbs <- update(snp1_model, n.iter=1000)
snp1_model_test <- coda.samples(snp1_model, c("beta_snp1", "OR1"), n.iter=10000)
summary(snp1_model_test)

#Gelman and Rubin Diagnostics
geweke.diag(snp1_model_test, frac1=0.1, frac2=0.5)
gelman.diag(snp1_model_test)
gelman.plot(snp1_model_test)
```



### SNP2 Analysis
```{r}
snp2_model_bugs <- 
"model{
#Logistic Model
for (i in 1:N){
  outcome[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_sex*sex[i] + beta_pc1*(pc1[i]-pc1_mean) +  beta_pc2*(pc2[i]-pc2_mean) +  beta_snp2*(snp2[i]) + beta_family[family[i]]
}

#Random Effect Per Family
for (i in 1:n_family){
  beta_family[i] ~ dnorm(beta0, tau)
}

#Prior Distributions
beta0 ~ dnorm(0, 0.01)
beta_sex ~ dnorm(0, 0.01)
beta_pc1 ~ dnorm(0, 0.01)
beta_pc2 ~ dnorm(0, 0.01)
beta_snp2 ~ dnorm(0, 0.01)
tau ~ dgamma(1, 1)

OR2 <- exp(beta_snp2)
pc1_mean <- mean(pc1[])
pc2_mean <- mean(pc2[])
}"

snp2_model <- jags.model(textConnection(snp2_model_bugs), data=SNP_data, n.adapt=1000, n.chains=3)
snp2_model_gibbs <- update(snp2_model, n.iter=1000)
snp2_model_test <- coda.samples(snp2_model, c("beta_snp2", "OR2"), n.iter=10000)
summary(snp2_model_test)

#Gelman and Rubin Diagnostics
geweke.diag(snp2_model_test, frac1=0.1, frac2=0.5)
gelman.diag(snp2_model_test)
gelman.plot(snp2_model_test)
```



### SNP3 Analysis
```{r}
snp3_model_bugs <- 
"model{
#Logistic Model
for (i in 1:N){
  outcome[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_sex*sex[i] + beta_pc1*(pc1[i]-pc1_mean) +  beta_pc2*(pc2[i]-pc2_mean) +  beta_snp3*(snp3[i]) + beta_family[family[i]]
}

#Random Effect Per Family
for (i in 1:n_family){
  beta_family[i] ~ dnorm(beta0, tau)
}

#Prior Distributions
beta0 ~ dnorm(0, 0.01)
beta_sex ~ dnorm(0, 0.01)
beta_pc1 ~ dnorm(0, 0.01)
beta_pc2 ~ dnorm(0, 0.01)
beta_snp3 ~ dnorm(0, 0.01)
tau ~ dgamma(1, 1)

OR3 <- exp(beta_snp3)
pc1_mean <- mean(pc1[])
pc2_mean <- mean(pc2[])
}"

snp3_model <- jags.model(textConnection(snp3_model_bugs), data=SNP_data, n.adapt=1000, n.chains=3)
snp3_model_gibbs <- update(snp3_model, n.iter=1000)
snp3_model_test <- coda.samples(snp3_model, c("beta_snp3", "OR3"), n.iter=10000)
summary(snp3_model_test)

#Gelman and Rubin Diagnostics
geweke.diag(snp3_model_test, frac1=0.1, frac2=0.5)
gelman.diag(snp3_model_test)
gelman.plot(snp3_model_test)
```



### SNP4 Analysis
```{r}
snp4_model_bugs <- 
"model{
#Logistic Model
for (i in 1:N){
  outcome[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_sex*sex[i] + beta_pc1*(pc1[i]-pc1_mean) +  beta_pc2*(pc2[i]-pc2_mean) +  beta_snp4*(snp4[i]) + beta_family[family[i]]
}

#Random Effect Per Family
for (i in 1:n_family){
  beta_family[i] ~ dnorm(beta0, tau)
}

#Prior Distributions
beta0 ~ dnorm(0, 0.01)
beta_sex ~ dnorm(0, 0.01)
beta_pc1 ~ dnorm(0, 0.01)
beta_pc2 ~ dnorm(0, 0.01)
beta_snp4 ~ dnorm(0, 0.01)
tau ~ dgamma(1, 1)

OR4 <- exp(beta_snp4)
pc1_mean <- mean(pc1[])
pc2_mean <- mean(pc2[])
}"

snp4_model <- jags.model(textConnection(snp4_model_bugs), data=SNP_data, n.adapt=1000, n.chains=3)
snp4_model_gibbs <- update(snp4_model, n.iter=1000)
snp4_model_test <- coda.samples(snp4_model, c("beta_snp4", "OR4"), n.iter=10000)
summary(snp4_model_test)

#Gelman and Rubin Diagnostics
geweke.diag(snp4_model_test, frac1=0.1, frac2=0.5)
gelman.diag(snp4_model_test)
gelman.plot(snp4_model_test)
```


### SNP5 Analysis
```{r}
snp5_model_bugs <- 
"model{
#Logistic Model
for (i in 1:N){
  outcome[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_sex*sex[i] + beta_pc1*(pc1[i]-pc1_mean) +  beta_pc2*(pc2[i]-pc2_mean) +  beta_snp5*(snp5[i]) + beta_family[family[i]]
}

#Random Effect Per Family
for (i in 1:n_family){
  beta_family[i] ~ dnorm(beta0, tau)
}

#Prior Distributions
beta0 ~ dnorm(0, 0.01)
beta_sex ~ dnorm(0, 0.01)
beta_pc1 ~ dnorm(0, 0.01)
beta_pc2 ~ dnorm(0, 0.01)
beta_snp5 ~ dnorm(0, 0.01)
tau ~ dgamma(1, 1)

OR5 <- exp(beta_snp5)
pc1_mean <- mean(pc1[])
pc2_mean <- mean(pc2[])
}"

snp5_model <- jags.model(textConnection(snp5_model_bugs), data=SNP_data, n.adapt=1000, n.chains=3)
snp5_model_gibbs <- update(snp5_model, n.iter=1000)
snp5_model_test <- coda.samples(snp5_model, c("beta_snp5", "OR5"), n.iter=10000)
summary(snp5_model_test)

#Gelman and Rubin Diagnostics
geweke.diag(snp5_model_test, frac1=0.1, frac2=0.5)
gelman.diag(snp5_model_test)
gelman.plot(snp5_model_test)
```



### SNP6 Analysis
```{r}
snp6_model_bugs <- 
"model{
#Logistic Model
for (i in 1:N){
  outcome[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_sex*sex[i] + beta_pc1*(pc1[i]-pc1_mean) +  beta_pc2*(pc2[i]-pc2_mean) +  beta_snp6*(snp6[i]) + beta_family[family[i]]
}

#Random Effect Per Family
for (i in 1:n_family){
  beta_family[i] ~ dnorm(beta0, tau)
}

#Prior Distributions
beta0 ~ dnorm(0, 0.01)
beta_sex ~ dnorm(0, 0.01)
beta_pc1 ~ dnorm(0, 0.01)
beta_pc2 ~ dnorm(0, 0.01)
beta_snp6 ~ dnorm(0, 0.01)
tau ~ dgamma(1, 1)

OR6 <- exp(beta_snp6)
pc1_mean <- mean(pc1[])
pc2_mean <- mean(pc2[])
}"

snp6_model <- jags.model(textConnection(snp6_model_bugs), data=SNP_data, n.adapt=1000, n.chains=3)
snp6_model_gibbs <- update(snp6_model, n.iter=1000)
snp6_model_test <- coda.samples(snp6_model, c("beta_snp6", "OR6"), n.iter=10000)
summary(snp6_model_test)

#Gelman and Rubin Diagnostics
geweke.diag(snp6_model_test, frac1=0.1, frac2=0.5)
gelman.diag(snp6_model_test)
gelman.plot(snp6_model_test)
```



### SNP7 Analysis
```{r}
snp7_model_bugs <- 
"model{
#Logistic Model
for (i in 1:N){
  outcome[i] ~ dbin(p[i], 1)
  logit(p[i]) <- beta_sex*sex[i] + beta_pc1*(pc1[i]-pc1_mean) +  beta_pc2*(pc2[i]-pc2_mean) +  beta_snp7*(snp7[i]) + beta_family[family[i]]
}

#Random Effect Per Family
for (i in 1:n_family){
  beta_family[i] ~ dnorm(beta0, tau)
}

#Prior Distributions
beta0 ~ dnorm(0, 0.01)
beta_sex ~ dnorm(0, 0.01)
beta_pc1 ~ dnorm(0, 0.01)
beta_pc2 ~ dnorm(0, 0.01)
beta_snp7 ~ dnorm(0, 0.01)
tau ~ dgamma(1, 1)

OR7 <- exp(beta_snp7)
pc1_mean <- mean(pc1[])
pc2_mean <- mean(pc2[])
}"

snp7_model <- jags.model(textConnection(snp7_model_bugs), data=SNP_data, n.adapt=1000, n.chains=3)
snp7_model_gibbs <- update(snp7_model, n.iter=1000)
snp7_model_test <- coda.samples(snp7_model, c("beta_snp7", "OR7"), n.iter=10000)
summary(snp7_model_test)

#Gelman and Rubin Diagnostics
geweke.diag(snp7_model_test, frac1=0.1, frac2=0.5)
gelman.diag(snp7_model_test)
gelman.plot(snp7_model_test)
```



### Boxplots of OR 
```{r}
snp1_model_OR <- as.matrix(snp1_model_test[,1])
snp2_model_OR <- as.matrix(snp2_model_test[,1])
snp3_model_OR <- as.matrix(snp4_model_test[,1])
snp4_model_OR <- as.matrix(snp4_model_test[,1])
snp5_model_OR <- as.matrix(snp5_model_test[,1])
snp6_model_OR <- as.matrix(snp6_model_test[,1])
snp7_model_OR <- as.matrix(snp7_model_test[,1])

snp_OR <- cbind(snp1_model_OR, snp2_model_OR, snp3_model_OR, snp4_model_OR, snp5_model_OR, snp6_model_OR, snp7_model_OR)
boxplot(snp_OR, col="hotpink", main="Boxplot of SNP ORs for Longevity")
```























