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
