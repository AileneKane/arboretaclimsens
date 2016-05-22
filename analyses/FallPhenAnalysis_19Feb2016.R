##Getting 2014 Fall and 2015 Spring phenology data in order
####Getting peak leaf color and peak leaf fall from my data
setwd("/Users/aileneettinger/Dropbox/Documents/Work/Putnam/Data_Ailene_18Feb2016_phylo")
fallphendat=read.csv("PhenData_Fall2014_10Dec.csv")
head(fallphendat)
colnames(fallphendat)
fallphendat=fallphendat[1:201,]
dim(fallphendat)
#rehapedata

fallphen.long = reshape(fallphendat, varying = list(names(fallphendat)[4:12], names(fallphendat)[13:21],names(fallphendat)[31:39],names(fallphendat)[40:48],names(fallphendat)[49:57]), direction = "long", v.names = c("Col_lvs_perc","Col_lvs", "Date","Fall_lvs_perc", "Fall_lvs"), times = c(5:13))
head(fallphen.long)
fallphen<-subset(fallphen.long,select=c("ACC_NUM","NAME","Date","Col_lvs","Col_lvs_perc", "Fall_lvs","Fall_lvs_perc"))
fallphen$Date<-as.Date(fallphen$Date,"%m/%d/%Y")
inds<-unique(fallphen$ACC_NUM)
firstcolall=c(rep(NA,times=length(inds)))
peakcolall=c(rep(NA,times=length(inds)))
species=c(rep(NA,times=length(inds)))
for(i in 1:length(inds)){
  inddat=fallphen[fallphen$ACC_NUM==inds[i],]
  peakcol1=inddat[inddat$Col_lvs_perc==max(inddat$Col_lvs_perc, na.rm = TRUE),]$Date
  if(is.na(peakcol1)[1]==TRUE) peakcol2=as.character(peakcol1[-which(is.na(peakcol1))][1])
  if(is.na(peakcol1)[1]==FALSE) peakcol2=as.character(peakcol1[1])
  firstcol1=min(inddat[inddat$Col_lvs=="Y",]$Date, na.rm = TRUE)
  species[i]=as.character(inddat$NAME[1])
  peakcolall[i]=peakcol2
  firstcolall[i]=as.character(firstcol1)
  }
peakfallall=c(rep(NA,times=length(inds)))
firstfallall=c(rep(NA,times=length(inds)))
for(i in 1:length(inds)){
  inddat=fallphen[fallphen$ACC_NUM==inds[i],]
  peakfall1=inddat[inddat$Fall_lvs_perc==max(inddat$Fall_lvs_perc, na.rm = TRUE),]$Date
  percfall1=inddat[inddat$Fall_lvs_perc==max(inddat$Fall_lvs_perc, na.rm = TRUE),]$Fall_lvs_perc
  if(is.na(peakfall1)[1]==TRUE) peakfall2=as.character(peakfall1[-which(is.na(peakfall1))][1])
  if(is.na(peakfall1)[1]==FALSE) peakfall2=as.character(peakfall1[1])
  if(is.na(percfall1)[1]==TRUE) percfall2=as.character(percfall1[-which(is.na(percfall1))][1])
  if(is.na(percfall1)[1]==FALSE) percfall2=as.character(percfall1[1])
  firstfall=min(inddat[inddat$Fall_lvs=="Y",]$Date, na.rm = TRUE)
  peakfallall[i]=peakfall2
  percfallall[i]=percfall2
  firstfallall[i]=as.character(firstfall)
}
allfallphendat=cbind(as.character(inds),species,firstcolall,peakcolall,firstfallall,peakfallall,percfallall)
colnames(allfallphendat)=c("Accession","Species","FirstColorDate","PeakColorDate","FirstDropDate","PeakDropDate","PercDropped")
write.csv(allfallphendat,"FallPhenDates.csv")
head(allfallphendat)

fallphendat<-read.csv("FallPhenDates.csv", header=T)
head(fallphendat)
fallphendat$PeakColDOY<- as.POSIXlt(fallphendat$PeakColorDate)$yday
fallphendat$FirstColDOY<- as.POSIXlt(fallphendat$FirstColorDate)$yday
fallphendat$FirstDropDOY<- as.POSIXlt(fallphendat$FirstDropDate)$yday
fallphendat$PeakDropDOY<- as.POSIXlt(fallphendat$PeakDropDate)$yday
head(fallphendat)
write.csv(fallphendat,"FallPhenDatesDoy.csv")

##Getting leaf budburst, leaf out, from my 2015 spring data
springphendat=read.csv("SpringPhenology2015.csv", header=T)
head(springphendat)
colnames(springphendat)
dim(springphendat)
springphendat<-springphendat[1:341,]
#rehapedata so that there are multiple rows for each individual (with a row for each survey date)
springphen.long = reshape(springphendat, varying = list(names(springphendat)[6:16], names(springphendat)[17:27],names(springphendat)[28:38],names(springphendat)[39:49],names(springphendat)[50:60]), direction = "long", v.names = c("BLB_NE","Date","FullExpLvs","ILS_NA","L_YN"))
head(springphen.long)
#convert data for budburst (BLB_NE), and leaf-out (leaves >95% in size), into all N and Y.
unique(springphen.long$BLB_NE)
springphen.long$Budburst=NA
springphen.long[which(springphen.long$BLB_NE==">10000"|springphen.long$BLB_NE=="1001-10000"|springphen.long$BLB_NE=="101-1000"|springphen.long$BLB_NE=="11-100"|springphen.long$BLB_NE=="3-10"|springphen.long$BLB_NE=="10-Mar"),]$Budburst<-"Y"
springphen.long[which(springphen.long$BLB_NE=="<3"),]$Budburst="N"
unique(springphen.long$Budburst)

springphen.long$LeafOut=NA
unique(springphen.long$FullExpLvs)#leafout based on when there is one fully expanded leaf noted
springphen.long[which(springphen.long$FullExpLvs==">10000"|springphen.long$FullExpLvs=="11-100"|springphen.long$FullExpLvs=="101-1000"|springphen.long$FullExpLvs=="1001-10000"|springphen.long$FullExpLvs=="3-10"|springphen.long$FullExpLvs=="3-11"),]$LeafOut="Y"
springphen.long[which(springphen.long$FullExpLvs=="<3"),]$LeafOut="N"
unique(springphen.long$LeafOut)
springphen.long[which(is.na(springphen.long$LeafOut)),]

springphen.long$LeafEmerg=NA
unique(springphen.long$L_YN)#number of leaves that have emerged
springphen.long[which(springphen.long$L_YN==">10000"|springphen.long$L_YN=="75-94"|springphen.long$L_YN=="101-1000"|springphen.long$L_YN=="1001-10000"|springphen.long$L_YN=="10-Mar"|springphen.long$L_YN=="3-10"|springphen.long$L_YN=="74-94"|springphen.long$L_YN=="25-49"|springphen.long$L_YN=="24-May"|springphen.long$L_YN=="50-74"|springphen.long$L_YN=="5 - 24"|springphen.long$L_YN==">95"|springphen.long$L_YN=="75-95"),]$LeafEmerg="Y"
springphen.long[which(springphen.long$FullExpLvs=="<3"),]$LeafOut="N"
unique(springphen.long$LeafEmerg)
springphen.long[which(is.na(springphen.long$LeafEmerg)),]

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
