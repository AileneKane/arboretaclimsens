#look at phylogenetic relationships in traits and climate sensitivity of arnold arboreutm tree cores!
#started december 7, 2016
#by Ailene Ettinger
rm(list=ls())
options(stringsAsFactors=FALSE)

library(dplR)
library(dplyr)
library(plyr)
library(lattice)
library(lme4)
library(RColorBrewer)
library(ape)
library(brms)
library(MCMCglmm)
#update.packages()
#load trait files
setwd("~/git/arboretaclimsens")
la <- read.csv("data/AA2015_leaf_areas.csv")
lm<-read.csv("data/AA2015_leaf_mass.csv")
woodmass<-read.csv("data/AA2015_twig_mass.csv")
height<-read.csv("data/TraitCollectionDatasheet2015_cleaned.csv", header=T)
springphen<-read.csv("data/SpringPhenDatesDoy.csv", header=T)
##add phenology data
fallphen<-read.csv("data/FallPhenDatesDoy.csv", header=T)#need to fix the names of fall phen to match phylogeny

#combine and clean trait files
source("cleaning/traitcleaning.R")
head(traits)
traits2<-na.omit(traits)#only 211 individuals have everything- need to clean/impute more data
#Add columns for the species mean of each trait
traits2$lma_sp<-with(traits2, sapply(split(lma, sp), mean)[sp])
traits2$dens_sp<-with(traits2, sapply(split(dens, sp), mean)[sp])
traits2$ht_sp<-with(traits2, sapply(split(ht, sp), mean)[sp])
traits2$phylo<-traits2$sp
head(traits2)
#read in phylogeny
aatree <- read.nexus("data/Ailene_comb.tre")
# Remove species in the tree that are not in the trait matrix
species.to.exclude <- aatree$tip.label[!(aatree$tip.label %in% 
                                           traits2$sp)]
phylo <- drop.tip(aatree,species.to.exclude)
rm(species.to.exclude)
#Remove species not in the phylogeny
traits3<-traits2[traits2$sp %in% phylo$tip.label,]

#prep the phylogeny 
inv.phylo <- MCMCglmm::inverseA(phylo, nodes = "TIPS", scale = TRUE)
A <- solve(inv.phylo$Ainv)
rownames(A) <- rownames(inv.phylo$Ainv)

#We want to model phenology as a function of traits (height, wood density, lma)
#budburst vs lma
bb_lma.mod<- brm(bb ~ lma_sp + (1|phylo) + (1|sp), 
                     data = traits3, family = gaussian(), 
                     cov_ranef = list(phylo = A),
                    sample_prior = TRUE, chains = 2, cores = 2, 
                     iter = 4000, warmup = 1000)

summary(bb_lma.mod)## + effect of lma, 356 divergent transitions
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(bb_lma.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(bb_lma.mod, hyp, class = NULL))
plot(hyp)
#update model with variability
traits3$lma_spvar <- traits3$lma - traits3$lma_sp
bb_lma.mod2 <- update(bb_lma.mod, formula = ~ . + lma_spvar,
                        newdata = traits3, chains = 2, cores = 2, 
                        iter = 4000, warmup = 1000)
summary(bb_lma.mod2)# 
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(bb_lma.mod2), points = TRUE) 
bb_lma.mod3<- brm(bb ~ lma + (1|phylo) + (1|sp), 
                 data = traits3, family = gaussian(), 
                 cov_ranef = list(phylo = A),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)
summary(bb_lma.mod3)#no divergent transitions, no effect of lma
plot(marginal_effects(bb_lma.mod), points = TRUE, bg=as.numeric(as.factor(traits3$group)))

#budburst vs wood dens
bb_dens.mod3<- brm(bb ~ dens + (1|phylo) + (1|sp), 
                 data = traits3, family = gaussian(), 
                 cov_ranef = list(phylo = A),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(bb_dens.mod3)#positive effect of density
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(bb_dens.mod3), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(bb_dens.mod3, hyp, class = NULL))
plot(hyp)# 0.89 
#update with within species variability
traits3$dens_spvar <- traits3$dens - traits3$dens_sp
bb_dens.mod2 <- update(bb_dens.mod, formula = ~ . + dens_spvar,
                      newdata = traits3, chains = 2, cores = 2, 
                      iter = 4000, warmup = 1000)
summary(bb_dens.mod2)#
plot(marginal_effects(bb_dens.mod2), points = TRUE) 
cols=c("lightblue", "darkgreen")
plot(traits3$dens,traits3$bb,pch=21, bg=cols[as.numeric(as.factor(traits3$group))])
plot(traits3$lma,traits3$bb,pch=21, bg=cols[as.numeric(as.factor(traits3$group))])
plot(traits3$ht,traits3$bb,pch=21, bg=cols[as.numeric(as.factor(traits3$group))])

#budburst vs height
bb_ht.mod3<-brm(bb ~ ht + (1|phylo) + (1|sp), 
                  data = traits3, family = gaussian(), 
                  cov_ranef = list(phylo = A),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(bb_ht.mod3)#no effect of
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(bb_ht.mod3), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(bb_ht.mod3, hyp, class = NULL))
plot(hyp)
#update with species variability
traits3$ht_spvar <- traits3$ht - traits3$ht_sp
bb_ht.mod2 <- update(bb_ht.mod, formula = ~ . + ht_spvar,
                       newdata = traits3, chains = 2, cores = 2, 
                       iter = 4000, warmup = 1000)
summary(bb_ht.mod2)#
plot(marginal_effects(bb_dens.mod2), points = TRUE) 


###now leaf emergence

le_lma.mod<- brm(le ~ lma_sp + (1|phylo) + (1|sp), 
                 data = traits3, family = gaussian(), 
                 cov_ranef = list(phylo = A),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(le_lma.mod)#lma has a positive effect
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(le_lma.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(le_lma.mod, hyp, class = NULL))
plot(hyp)
#update model with variability
traits3$lma_spvar <- traits3$lma - traits3$lma_sp
le_lma.mod2 <- update(le_lma.mod, formula = ~ . + lma_spvar,
                      newdata = traits3, chains = 2, cores = 2, 
                      iter = 4000, warmup = 1000)
summary(le_lma.mod2)#
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(le_lma.mod2), points = TRUE) 

le_lma.mod3<- brm(le ~ lma + (1|phylo) + (1|sp), 
                  data = traits3, family = gaussian(), 
                  cov_ranef = list(phylo = A),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(le_lma.mod3)#lma has a positive effect
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(le_lma.mod3), points = TRUE) 

#leaf emergence vs wood dens
le_dens.mod3<- brm(le ~ dens+ (1|phylo) + (1|sp), 
                  data = traits3, family = gaussian(), 
                  cov_ranef = list(phylo = A),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(le_dens.mod3)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(le_dens.mod3), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(le_dens.mod, hyp, class = NULL))
plot(hyp)

#update model with variability
le_dens.mod2 <- update(le_dens.mod, formula = ~ . + dens_spvar,
                      newdata = traits3, chains = 2, cores = 2, 
                      iter = 4000, warmup = 1000)
summary(le_dens.mod2)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(le_dens.mod2), points = TRUE) 


#leaf emergence vs height
le_ht.mod3<- brm(le ~ ht + (1|phylo) + (1|sp), 
                data = traits3, family = gaussian(), 
                cov_ranef = list(phylo = A),
                sample_prior = TRUE, chains = 2, cores = 2, 
                iter = 4000, warmup = 1000)

summary(le_ht.mod3)#
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(le_ht.mod3), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(le_ht.mod, hyp, class = NULL))
plot(hyp)


#update model with variability
le_ht.mod2 <- update(le_ht.mod, formula = ~ . + ht_spvar,
                      newdata = traits3, chains = 2, cores = 2, 
                      iter = 4000, warmup = 1000)
summary(le_ht.mod2)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(le_ht.mod2), points = TRUE) 
#Fit a model with all three traits
bb_all.mod3<-brm(bb ~ lma+dens+ht + (1|phylo) + (1|sp), 
                data = traits3, family = gaussian(), 
                cov_ranef = list(phylo = A),
                sample_prior = TRUE, chains = 2, cores = 2, 
                iter = 4000, warmup = 1000)
summary(bb_all.mod3)#
plot(marginal_effects(bb_all.mod3), points = TRUE) 
le_all.mod3<-brm(le ~ lma+dens+ht + (1|phylo) + (1|sp), 
                 data = traits3, family = gaussian(), 
                 cov_ranef = list(phylo = A),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(le_all.mod3)#
plot(marginal_effects(le_all.mod3), points = TRUE) 


#Now analyze gymnosperms and angiosperms separately
gtraits<-traits3[traits3$group=="gymnosperm",]
atraits<-traits3[traits3$group=="angiosperm",]
gspecies.to.exclude <- aatree$tip.label[!(aatree$tip.label %in% 
                                           gtraits$sp)]
aspecies.to.exclude <- aatree$tip.label[!(aatree$tip.label %in% 
                                            atraits$sp)]

gphylo <- drop.tip(aatree,gspecies.to.exclude)
aphylo <- drop.tip(aatree,aspecies.to.exclude)
gtraits2<-gtraits[gtraits$sp %in% gphylo$tip.label,]
atraits2<-atraits[atraits$sp %in% aphylo$tip.label,]

#prep the gymnosperm phylogeny 
inv.gphylo <- MCMCglmm::inverseA(gphylo, nodes = "TIPS", scale = TRUE)
gA <- solve(inv.gphylo$Ainv)
rownames(gA) <- rownames(inv.gphylo$Ainv)

#prep the angiosperm phylogeny 
inv.aphylo <- MCMCglmm::inverseA(aphylo, nodes = "TIPS", scale = TRUE)
aA <- solve(inv.aphylo$Ainv)
rownames(aA) <- rownames(inv.aphylo$Ainv)

#bb
gbb_lma.mod3<- brm(bb ~ lma + (1|phylo) + (1|sp), 
                 data = gtraits2, family = gaussian(), 
                 cov_ranef = list(phylo = gA),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(gbb_lma.mod3)## + effect of lma, 356 divergent transitions
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(gbb_lma.mod3), points = TRUE) 
#gtraits$lma_spvar <- gtraits$lma - gtraits$lma_sp

#gbb_lma.mod2 <- update(gbb_lma.mod, formula = ~ . + lma_spvar,
 #                      newdata = gtraits, chains = 2, cores = 2, 
    #                   iter = 4000, warmup = 1000)
#summary(gbb_lma.mod2)
#plot(marginal_effects(gbb_lma.mod2), points = TRUE) 


abb_lma.mod<- brm(bb ~ lma + (1|phylo) + (1|sp), 
                  data = atraits2, family = gaussian(), 
                  cov_ranef = list(phylo = aA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(abb_lma.mod)## + effect of lma, 356 divergent transitions
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(abb_lma.mod), points = TRUE) 

#atraits$lma_spvar <- atraits$lma - atraits$lma_sp
#abb_lma.mod2 <- update(abb_lma.mod, formula = ~ . + lma_spvar,
#                       newdata = atraits, chains = 2, cores = 2, 
                       iter = 4000, warmup = 1000)


#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(abb_lma.mod, hyp, class = NULL))
plot(hyp)
#update model with variability
#traits3$lma_spvar <- traits3$lma - traits3$lma_sp
#bb_lma.mod2 <- update(bb_lma.mod, formula = ~ . + lma_spvar,
                      newdata = traits3, chains = 2, cores = 2, 
                      iter = 4000, warmup = 1000)
#summary(bb_lma.mod)# + effect of lma, 946 divergent transitions
#plot(bbmod, N = 2, ask = FALSE)
#plot(marginal_effects(bb_lma.mod2), points = TRUE) 


#budburst vs wood dens
gbb_dens.mod<- brm(bb ~ dens + (1|phylo) + (1|sp), 
                  data = gtraits2, family = gaussian(), 
                  cov_ranef = list(phylo = gA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(gbb_dens.mod)#positive effect of density
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(gbb_dens.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gbb_dens.mod, hyp, class = NULL))
plot(hyp)
#update with within species variability
#traits3$dens_spvar <- traits3$dens - traits3$dens_sp
#bb_dens.mod2 <- update(bb_dens.mod, formula = ~ . + dens_spvar,
#                       newdata = traits3, chains = 2, cores = 2, 
 #                      iter = 4000, warmup = 1000)
#summary(bb_dens.mod2)#
#plot(marginal_effects(bb_dens.mod2), points = TRUE) 

#budburst vs height
gbb_ht.mod<- brm(bb ~ ht + (1|phylo) + (1|sp), 
                data = gtraits2, family = gaussian(), 
                cov_ranef = list(phylo = gA),
                sample_prior = TRUE, chains = 2, cores = 2, 
                iter = 4000, warmup = 1000)

summary(gbb_ht.mod)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(gbb_ht.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gbb_ht.mod, hyp, class = NULL))
plot(hyp)
#update with species variability
#traits3$ht_spvar <- traits3$ht - traits3$ht_sp
#bb_ht.mod2 <- update(bb_ht.mod, formula = ~ . + ht_spvar,
                     newdata = traits3, chains = 2, cores = 2, 
                     iter = 4000, warmup = 1000)
#summary(bb_ht.mod2)#
#plot(marginal_effects(bb_dens.mod2), points = TRUE) 


###now leaf emergence
gle_lma.mod<- brm(le ~ lma + (1|phylo) + (1|sp), 
                 data = gtraits2, family = gaussian(), 
                 cov_ranef = list(phylo = gA),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(gle_lma.mod)#lma has a positive effect
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(gle_lma.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gle_lma.mod, hyp, class = NULL))
plot(hyp)
#update model with variability
#le_lma.mod2 <- update(le_lma.mod, formula = ~ . + lma_spvar,
       #               newdata = traits3, chains = 2, cores = 2, 
        #              iter = 4000, warmup = 1000)
#summary(le_lma.mod2)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
#plot(marginal_effects(le_lma.mod2), points = TRUE) 


#leaf emergence vs wood dens
ale_dens.mod<- brm(le ~ dens + (1|phylo) + (1|sp), 
                  data = atraits2, family = gaussian(), 
                  cov_ranef = list(phylo = aA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(ale_dens.mod)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(ale_dens.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(ale_dens.mod, hyp, class = NULL))
plot(hyp)

#update model with variability
#le_dens.mod2 <- update(le_dens.mod, formula = ~ . + dens_spvar,
 #                      newdata = gtraits2, chains = 2, cores = 2, 
  #                     iter = 4000, warmup = 1000)
#summary(le_dens.mod2)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
#plot(marginal_effects(le_dens.mod2), points = TRUE) 


#leaf emergence vs height
gle_ht.mod<- brm(le ~ ht + (1|phylo) + (1|sp), 
                data = gtraits2, family = gaussian(), 
                cov_ranef = list(phylo = gA),
                sample_prior = TRUE, chains = 2, cores = 2, 
                iter = 4000, warmup = 1000)

summary(gle_ht.mod)#
#plot(bbmod, N = 2, ask = FALSE)
plot(marginal_effects(gle_ht.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gle_ht.mod, hyp, class = NULL))
plot(hyp)


#update model with variability
#le_ht.mod2 <- update(le_ht.mod, formula = ~ . + ht_spvar,
#                     newdata = traits3, chains = 2, cores = 2, 
#                     iter = 4000, warmup = 1000)
#summary(le_ht.mod2)#no effect of lma
#plot(bbmod, N = 2, ask = FALSE)
#plot(marginal_effects(le_ht.mod2), points = TRUE) 


#prep the phylogenies 
gspecies.to.exclude <- aatree$tip.label[!(aatree$tip.label %in% 
                                           gtraits$sp)]
gphylo <- drop.tip(aatree,gspecies.to.exclude)
aspecies.to.exclude <- aatree$tip.label[!(aatree$tip.label %in% 
                                            atraits$sp)]
aphylo <- drop.tip(aatree,aspecies.to.exclude)

ginv.phylo <- MCMCglmm::inverseA(gphylo, nodes = "TIPS", scale = TRUE)
ainv.phylo <- MCMCglmm::inverseA(aphylo, nodes = "TIPS", scale = TRUE)

gA <- solve(ginv.phylo$Ainv)
rownames(gA) <- rownames(ginv.phylo$Ainv)

aA <- solve(ainv.phylo$Ainv)
rownames(aA) <- rownames(ainv.phylo$Ainv)

#gymnosperm budburst vs lma
gbb_lma.mod<- brm(bb ~ lma + (1|phylo) + (1|sp), 
                 data = gtraits, family = gaussian(), 
                 cov_ranef = list(phylo = gA),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(gbb_lma.mod)## + effect of lma, 648 divergent transitions
plot(marginal_effects(gbb_lma.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gbb_lma.mod, hyp, class = NULL))
plot(hyp)
#angiosperm budburst vs lma
abb_lma.mod<- brm(bb ~ lma + (1|phylo) + (1|sp), 
                  data = atraits, family = gaussian(), 
                  cov_ranef = list(phylo = aA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(abb_lma.mod)## + effect of lma, 356 divergent transitions
plot(marginal_effects(abb_lma.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(abb_lma.mod, hyp, class = NULL))
plot(hyp)

#gymnosperm budburst vs dens
gbb_dens.mod<- brm(bb ~ dens + (1|phylo) + (1|sp), 
                  data = gtraits, family = gaussian(), 
                  cov_ranef = list(phylo = gA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(gbb_dens.mod)## + effect of lma, 648 divergent transitions
plot(marginal_effects(gbb_dens.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gbb_dens.mod, hyp, class = NULL))
plot(hyp)
#angiosperm budburst vs dens
abb_dens.mod<- brm(bb ~ dens + (1|phylo) + (1|sp), 
                  data = atraits, family = gaussian(), 
                  cov_ranef = list(phylo = aA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(abb_dens.mod)## + effect of lma, 356 divergent transitions
plot(marginal_effects(abb_dens.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(abb_dens.mod, hyp, class = NULL))
plot(hyp)


#gymnosperm budburst vs lma
gbb_ht.mod<- brm(bb ~ ht + (1|phylo) + (1|sp), 
                   data = gtraits, family = gaussian(), 
                   cov_ranef = list(phylo = gA),
                   sample_prior = TRUE, chains = 2, cores = 2, 
                   iter = 4000, warmup = 1000)

summary(gbb_ht.mod)## + effect of lma, 648 divergent transitions
plot(marginal_effects(gbb_ht.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gbb_ht.mod, hyp, class = NULL))
plot(hyp)
#angiosperm budburst vs lma
abb_ht.mod<- brm(bb ~ ht + (1|phylo) + (1|sp), 
                   data = atraits, family = gaussian(), 
                   cov_ranef = list(phylo = aA),
                   sample_prior = TRUE, chains = 2, cores = 2, 
                   iter = 4000, warmup = 1000)

summary(abb_ht.mod)## + effect of lma, 356 divergent transitions
plot(marginal_effects(abb_ht.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(abb_ht.mod, hyp, class = NULL))
plot(hyp)

#now leaf emergence

#gymnosperm le vs lma
gle_lma.mod<- brm(le ~ lma+ (1|phylo) + (1|sp), 
                  data = gtraits2, family = gaussian(), 
                  cov_ranef = list(phylo = gA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(gle_lma.mod)## + effect of lma, 648 divergent transitions
plot(marginal_effects(gle_lma.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gle_lma.mod, hyp, class = NULL))
plot(hyp)
#angiosperm budburst vs lma
ale_lma.mod<- brm(le ~ lma + (1|phylo) + (1|sp), 
                  data = atraits2, family = gaussian(), 
                  cov_ranef = list(phylo = aA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)

summary(ale_lma.mod)## + effect of lma, 356 divergent transitions
plot(marginal_effects(ale_lma.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(ale_lma.mod, hyp, class = NULL))
plot(hyp)

#gymnosperm budburst vs dens
gle_dens.mod<- brm(le ~ dens + (1|phylo) + (1|sp), 
                   data = gtraits2, family = gaussian(), 
                   cov_ranef = list(phylo = gA),
                   sample_prior = TRUE, chains = 2, cores = 2, 
                   iter = 4000, warmup = 1000)

summary(gle_dens.mod)## + effect of lma, 648 divergent transitions
plot(marginal_effects(gle_dens.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gle_dens.mod, hyp, class = NULL))
plot(hyp)
#angiosperm budburst vs dens
ale_dens.mod<- brm(le~ dens + (1|phylo) + (1|sp), 
                   data = atraits2, family = gaussian(), 
                   cov_ranef = list(phylo = aA),
                   sample_prior = TRUE, chains = 2, cores = 2, 
                   iter = 4000, warmup = 1000)

summary(ale_dens.mod)## + effect of lma, 356 divergent transitions
plot(marginal_effects(ale_dens.mod), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(ale_dens.mod, hyp, class = NULL))
plot(hyp)


#gymnosperm budburst vs lma
gle_ht.mod<- brm(le ~ ht + (1|phylo) + (1|sp), 
                 data = gtraits2, family = gaussian(), 
                 cov_ranef = list(phylo = gA),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(gle_ht.mod)## + effect of lma, 648 divergent transitions
plot(marginal_effects(gle_ht.mod), points = TRUE) 
#gle_ht.mod2 <- update(gle_ht.mod, formula = ~ . + ht_spvar,
 #                      newdata = gtraits, chains = 2, cores = 2, 
  #                     iter = 4000, warmup = 1000)
#summary(gle_ht.mod2)
#plot(marginal_effects(gle_ht.mod2), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(gle_ht.mod, hyp, class = NULL))
plot(hyp)
#angiosperm budburst vs lma
ale_ht.mod<- brm(le ~ ht + (1|phylo) + (1|sp), 
                 data = atraits2, family = gaussian(), 
                 cov_ranef = list(phylo = aA),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)

summary(ale_ht.mod)## + effect of lma, 356 divergent transitions
plot(marginal_effects(ale_ht.mod), points = TRUE) 
ale_ht.mod2 <- update(ale_ht.mod, formula = ~ . + ht_spvar,
                       newdata = atraits, chains = 2, cores = 2, 
                       iter = 4000, warmup = 1000)
summary(ale_ht.mod2)
plot(marginal_effects(gbb_lma.mod2), points = TRUE) 

#test for phylogenetic signal
hyp <- "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0"
(hyp <- hypothesis(ale_ht.mod, hyp, class = NULL))
plot(hyp)

#Does spring phenology predict fall phenology?
allphen<-na.omit(allphen)#200 individuals left after
#aphylo$tip.label[!(aphylo$tip.label %in%  allphen$sp)]#all species present in angiosperm phylog

#Remove species not in the phylogeny
allphen2<-allphen[allphen$sp %in% aphylo$tip.label,]

#prep the phylogeny 
ainv.phylo <- MCMCglmm::inverseA(aphylo, nodes = "TIPS", scale = TRUE)
aA <- solve(inv.phylo$Ainv)
rownames(aA) <- rownames(ainv.phylo$Ainv)
allphen2$phylo<-allphen2$sp
#fit bb models
fcol_bb.mod<- brm(col ~ bb + (1|phylo) + (1|sp), 
                 data = allphen2, family = gaussian(), 
                 cov_ranef = list(phylo = aA),
                 sample_prior = TRUE, chains = 2, cores = 2, 
                 iter = 4000, warmup = 1000)
summary(fcol_bb.mod)
plot(marginal_effects(fcol_bb.mod), points = TRUE) 

fdrop_bb.mod<- brm(drop ~ bb + (1|phylo) + (1|sp), 
                  data = allphen2, family = gaussian(), 
                  cov_ranef = list(phylo = aA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)
summary(fdrop_bb.mod)
plot(marginal_effects(fdrop_bb.mod), points = TRUE) 

pkdrop_bb.mod<- brm(pkdrop ~ bb + (1|phylo) + (1|sp), 
                   data = allphen2, family = gaussian(), 
                   cov_ranef = list(phylo = aA),
                   sample_prior = TRUE, chains = 2, cores = 2, 
                   iter = 4000, warmup = 1000)
summary(pkdrop_bb.mod)#get lots of divergent transitions
plot(marginal_effects(pkdrop_bb.mod), points = TRUE) 
#now leaf emergence
fcol_le.mod<- brm(col ~ le + (1|phylo) + (1|sp), 
                  data = allphen2, family = gaussian(), 
                  cov_ranef = list(phylo = aA),
                  sample_prior = TRUE, chains = 2, cores = 2, 
                  iter = 4000, warmup = 1000)
summary(fcol_le.mod)
plot(marginal_effects(fcol_le.mod), points = TRUE) 

fdrop_le.mod<- brm(drop ~ le + (1|phylo) + (1|sp), 
                   data = allphen2, family = gaussian(), 
                   cov_ranef = list(phylo = aA),
                   sample_prior = TRUE, chains = 2, cores = 2, 
                   iter = 4000, warmup = 1000)
summary(fdrop_le.mod)
plot(marginal_effects(fdrop_le.mod), points = TRUE) 

pkdrop_le.mod<- brm(pkdrop ~ le + (1|phylo) + (1|sp), 
                    data = allphen2, family = gaussian(), 
                    cov_ranef = list(phylo = aA),
                    sample_prior = TRUE, chains = 2, cores = 2, 
                    iter = 4000, warmup = 1000)
summary(pkdrop_le.mod)#24 divergent transitions
plot(marginal_effects(pkdrop_le.mod), points = TRUE) 

