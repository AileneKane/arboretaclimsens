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

