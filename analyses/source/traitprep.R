#traitprep for analysis
traits2<-na.omit(traits)#only 211 individuals have everything- need to impute more data
#Add columns for the species mean of each trait
traits2$lma_sp<-with(traits2, sapply(split(lma, sp), mean)[sp])
traits2$dens_sp<-with(traits2, sapply(split(dens, sp), mean)[sp])
traits2$ht_sp<-with(traits2, sapply(split(ht, sp), mean)[sp])
traits2$phylo<-traits2$sp