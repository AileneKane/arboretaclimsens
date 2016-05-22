##Getting leaf budburst, leaf out, flowering times from my (sparse) 2014 spring data
setwd("~/Dropbox/Documents/Work/Putnam/Data")
springphendat=read.csv("PhenData_Spring2014_coredtrees.csv", header=T)
head(springphendat)
colnames(springphendat)
dim(springphendat)
#first, check the data for mistakes
#unique(springphendat$Date2)#did this for all Date columns
#unique(springphendat$BLB_NE9)#did this for all BLB_NE9 columns
#unique(springphendat$L_YN3)#did this for all L_YN9 columns
#unique(springphendat$ILS_NA9)#did this for all ILS_NA9 columns
#unique(springphendat$FFB_PC1)#did this for all FFB_PC1columns
#unique(springphendat$OF_OPC9)#did this for all OF_OPC1 rc columns

#rehapedata so that there are multiple rows for each individual (with a row for each survey date)
springphen.long = reshape(springphendat, varying = list(names(springphendat)[6:14], names(springphendat)[15:23],names(springphendat)[24:32],names(springphendat)[33:41],names(springphendat)[42:50], names(springphendat)[51:59]),direction = "long", v.names = c("Date","BLB_NE","FFB_PC","ILS_NA","L_YN","OF_OPC"))
head(springphen.long)
#convert data for budburst (BLB_NE), open flowers, and leaf-out (leaves >95% in size), into all N and Y.
unique(springphen.long$BLB_NE)
#isoinds[which(isoinds$Stage2=="ArrSm"|isoinds$Stage2=="FlowSm"),]$Stage3="ArFlSm"
springphen.long$Budburst=NA
springphen.long[which(springphen.long$BLB_NE==">10000"|springphen.long$BLB_NE==">95%"|springphen.long$BLB_NE=="1001-10000"|springphen.long$BLB_NE=="101-1000"|springphen.long$BLB_NE=="11-100"|springphen.long$BLB_NE=="3-10"|springphen.long$BLB_NE=="5-24%"|springphen.long$BLB_NE=="25-49%"|springphen.long$BLB_NE=="75-95%"|springphen.long$BLB_NE=="50-74%"|springphen.long$BLB_NE=="Y"),]$Budburst="Y"
springphen.long[which(springphen.long$BLB_NE=="<3"|springphen.long$BLB_NE=="<5%"|springphen.long$BLB_NE=="N"),]$Budburst="N"
unique(springphen.long$Budburst)
springphen.long$OpenFlrs=NA
unique(springphen.long$OF_OPC)
springphen.long[which(springphen.long$OF_OPC==">10000"|springphen.long$OF_OPC==">95%"|springphen.long$OF_OPC=="1001-10000"|springphen.long$OF_OPC=="25-49%"|springphen.long$OF_OPC=="5-24%"|springphen.long$OF_OPC=="50-74%"|springphen.long$OF_OPC=="75-95%"|springphen.long$OF_OPC=="past"|springphen.long$OF_OPC=="Y"),]$OpenFlrs="Y"
springphen.long[which(springphen.long$OF_OPC=="<3"|springphen.long$OF_OPC=="<5%"|springphen.long$OF_OPC=="N"|springphen.long$OF_OPC=="0"),]$OpenFlrs="N"
unique(springphen.long$OpenFlrs)
springphen.long$LeafOut=NA
unique(springphen.long$ILS)
springphen.long[which(springphen.long$ILS==">95%"|springphen.long$ILS=="25-49%"|springphen.long$ILS=="5-24%"|springphen.long$ILS=="50-74%"|springphen.long$BLB_NE=="75-95%"),]$LeafOut="Y"
springphen.long[which(springphen.long$ILS=="<3"|springphen.long$ILS=="<5%"|springphen.long$ILS=="N"),]$LeafOut="N"
unique(springphen.long$LeafOut)

inds<-unique(springphen.long$Accession)
firstflow=c(rep(NA,times=length(inds)))
leafout=c(rep(NA,times=length(inds)))
leafbudburst=c(rep(NA,times=length(inds)))
species=c(rep(NA,times=length(inds)))
for(i in 1:length(inds)){
  inddat=springphen.long[springphen.long$Accession==inds[i],]
  firstflow[i]=as.character(inddat[which(inddat$OpenFlrs=="Y"),]$Date[1])
  leafout[i]=as.character(inddat[which(inddat$LeafOut=="Y"),]$Date[1])
  leafbudburst[i]=as.character(inddat[which(inddat$Budburst=="Y"),]$Date[1])
  species[i]=as.character(inddat$Species[1])  
}
alldat=cbind(as.character(inds),species,leafbudburst,leafout,firstflow)
colnames(alldat)=c("Accession","Species","LeafBudburst","LeafOut","FirstOpenFlowers")
head(alldat)
write.csv(alldat,"SpringPhenDates.csv")
