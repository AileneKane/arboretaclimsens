#Quantifying climate sensitivity

growth<-read.csv("data/ArnArboringwidths_processed_gymnosperms.csv")
head(growth)

#Crossdating for loop: this prints the stats for flagged cores and images for ALL cores (i couldn't fgure out how to only print for flagged cores)to help you diagnose cores that are poorly crossdated and should be checked in windendro
#for (i in 1:length(species)){
#  spdat<-coredata[coredata$Species==species[i],]
#  for (j in 1:length(stands)){
#    print(stands[j]);print(species[i])
#    spstdat<-spdat[spdat$Stand==stands[j],]
#    spstdat.rwl<-t(spstdat[,10:226])
#    rownames(spstdat.rwl)<-c(2007:1791)
#    colnames(spstdat.rwl)<-spstdat$Tag
#    spstdat.rwl.df<-as.data.frame(spstdat.rwl)
#    quartz()
#    spst.stats<-corr.rwl.seg(spstdat.rwl,seg.length=50,main=paste(stands[j],species[i]))#seg.length can be changed to analyze cores on a finer scale (right now set at 50 years)		
#    print(spst.stats$flags)
#    flagged<-spst.stats$flags
#    for (k in 1:dim(spstdat.rwl)[2]){
#      spstdat.rwl<-t(spstdat[,10:226])
#      rownames(spstdat.rwl)<-c(2007:1791)
#      colnames(spstdat.rwl)<-spstdat$Tag
#      spstdat.rwl.df<-as.data.frame(spstdat.rwl)
#      rownames(spstdat.rwl.df)<-c(2007:1791)
#     colnames(spstdat.rwl.df)<-spstdat$Tag
#      flagtree<-colnames(spstdat.rwl.df)[k]
#      flagged=spstdat.rwl.df[,k]
#      names(flagged)=rownames(spstdat.rwl.df)
#     spstdat.rwl.df[,k]=NULL
#quartz(width=12,height=5)
#seg.20<-corr.series.seg(rwl=spstdat.rwl.df,series=flagged,seg.length=20,bin.floor=0, main=paste(flagtree))#you may want to look at these line plots, but they use up too many windows if you are doing it for a whole stand
#      quartz(width=12,height=5)
#      ccf.20=ccf.series.rwl(rwl=spstdat.rwl.df,series=flagged,seg.length=20,bin.floor=0,main=paste(colnames(spstdat.rwl.df)[k]))
#    }}}

#detrend tree growth
#first get rid of 1st and last year, which are not reliable
rawdat<-growth[,c(1:2,9:159)] #only columns we need:Treename, Species,Years 2013-1861
#for each row, remove first ring (which is inaccurately measured)
get.first.ring <- function(x) max(which(!is.na(x))) 
f<-apply(rawdat,1,get.first.ring)
inds<-unique(rawdat$TreeName)
for(i in 1:length(inds)){
  rawdat[i,f[i]]<-NA
}
inddat<-rawdat[,1:2]
rngdat<-as.data.frame(rawdat[3:53])#select only most recent 50 years
rownames(rngdat)<-rawdat[,1]
rngdat2<-as.data.frame(t(rngdat))
rwidat<-detrend(rngdat2,method="Spline", f=0.5)
#get Warning message:In get(var, frame, inherits = FALSE) : internal error -3 in R_decompress1
dat.rwi<-cbind(inddat,t(rwidat))
dat.rwi.long<-reshape(dat.rwi,varying=list(names(dat.rwi[3:53])),direction = "long",v.names=c("rwi"),times=2013:1963)
rownames(dat.rwi.long)<-NULL
colnames(dat.rwi.long)<-c("treeid","species","year","rwi")
rwi<-dat.rwi.long[,-5]
#summarize climate data: just use temperature for now
climate<-read.csv("data/clim_bluehills.csv", na.strings = "-9999",header=TRUE)
bluehillsclim<-climate[which(climate$STATION_NAME=="BLUE HILL MA US"),]
temp<-subset(bluehillsclim,select=c(DATE,TMAX,TMIN))
temp$year<-substr(temp$DATE,1,4)
temp$month<-substr(temp$DATE,5,6)
temp$day<-substr(temp$DATE,7,8)
temp$TMEAN<-(temp$TMAX+temp$TMIN)/2#daily mean temp
temp$year<-as.numeric(temp$year)
temp$growyear<-temp$year#growing year=nov-oct
temp[which(temp$month=="11"|temp$month=="12"),]$growyear<-temp[which(temp$month=="11"|temp$month=="12"),]$year+1
#right now this is daily climate data. summarize to annual values. for now, just MAT
nov<-subset(temp,month=="11")
dec<-subset(temp,month=="12")
jan<-subset(temp,month=="01")
feb<-subset(temp,month=="02")
mar<-subset(temp,month=="03")
apr<-subset(temp,month=="04")
may<-subset(temp,month=="05")
jun<-subset(temp,month=="06")
jul<-subset(temp,month=="07")
aug<-subset(temp,month=="08")
sep<-subset(temp,month=="09")
winter<-rbind(dec,jan,feb)
dorm<-rbind(nov,dec,jan,feb,mar)
summer<-rbind(jun,jul,aug)
grow<-rbind(may,jun,jul,aug,sep)
mwint<-aggregate(winter$TMEAN,by=list(winter$growyear),mean,na.rm=T)
mdorm<-aggregate(dorm$TMEAN,by=list(dorm$growyear),mean,na.rm=T)
mdorm<-mdorm[-which(mdorm$Group.1==2017),]
mwint<-mwint[-which(mwint$Group.1==2017),]
msum<-aggregate(summer$TMEAN,by=list(summer$growyear),mean,na.rm=T)
mgrow<-aggregate(grow$TMEAN,by=list(grow$growyear),mean,na.rm=T)
mat<-aggregate(temp$TMEAN,by=list(temp$growyear),mean,na.rm=T)
mat<-mat[-which(mat$Group.1==2017),]
allclim<-cbind(mat,mwint[,2],mdorm[,2],msum[,2],mgrow[,2])
colnames(allclim)<-c("year","mat","wint","dormt","sumt","growt")

allclim<- allclim[order(allclim$year,decreasing = TRUE),]
allclim2<-allclim[which(allclim$year<2014),]
rwi2<-rwi[which(rwi$year>=1893),]
rwi2<-rwi2[order(rwi2$species),]
alldat<-join(rwi2,allclim2,by="year",type="left",match="all")
species<-unique(alldat$species)
#fit models to get climate sensitivity by species
alldat$treeid<-as.factor(alldat$treeid)
alldat$year<-as.factor(alldat$year)
modsums<-as.data.frame(matrix(data=NA,nrow=length(species),ncol=13))
colnames(modsums)<-c("species","null","mat_int","mat_B","win_int","win_B","dorm_int","dorm_B","sum_int","sum_B","grow_int","grow_B","bestmod")
for (i in 1:length(species)){
  spdat<-alldat[alldat$species==species[i],]
  #quartz()
  #xyplot(spdat$rwi ~ spdat$mat| spdat$treeid,groups = spdat$treeid,pch=16,type =c("p","r"),main=paste(species[i]),lwd=2)
  if(length(unique(spdat$treeid))>1)
  {
    null<-lmer(rwi~1+(1|treeid),REML=FALSE, data=spdat)
    matm<-lmer(rwi~mat+(1|treeid),REML=FALSE, data=spdat)
    winm<-lmer(rwi~wint+(1|treeid),REML=FALSE, data=spdat)
    dorm<-lmer(rwi~dormt+(1|treeid),REML=FALSE, data=spdat)
    sum<-lmer(rwi~sumt+(1|treeid),REML=FALSE, data=spdat)
    growm<-lmer(rwi~growt+(1|treeid),REML=FALSE, data=spdat)
    aic<-AIC(null,matm,winm,dorm,sum,growm)
    bestmod<-rownames(aic[which(aic$AIC==min(aic$AIC)),])
    modsums[i,]<-c(paste(species[i]),round(fixef(null),digits=3),round(fixef(matm),digits=3),round(fixef(winm),digits=3),round(fixef(dorm),digits=3),round(fixef(sum),digits=3),round(fixef(growm),digits=3),paste(bestmod))
  }
  else
  { 
    null<-lm(rwi~1,data=spdat)
    matm<-lm(rwi~mat,data=spdat)
    winm<-lm(rwi~wint,data=spdat)
    dorm<-lm(rwi~dormt,data=spdat)
    sum<-lm(rwi~sumt,data=spdat)
    growm<-lm(rwi~growt,data=spdat)
    aic<-AIC(null,matm,winm,dorm,sum,growm)
    bestmod<-rownames(aic[which(aic$AIC==min(aic$AIC)),])
    modsums[i,]<-c(paste(species[i]),round(coef(null),digits=3),round(coef(matm),digits=3),round(coef(winm),digits=3),round(coef(dorm),digits=3),round(coef(sum),digits=3),round(coef(growm),digits=3),paste(bestmod))
    print(paste(species[i]));print(bestmod)
  }
}
length(which(modsums$bestmod=="null"))#14
length(which(modsums$bestmod=="sum"))#9
length(which(modsums$bestmod=="winm"))#9
length(which(modsums$bestmod=="dorm"))#5
length(which(modsums$bestmod=="growm"))#2
modsums[which(modsums$bestmod=="dorm"),]$dorm_B#4/5 positive
modsums[which(modsums$bestmod=="winm"),]$win_B#6/9 negative
modsums[which(modsums$bestmod=="sum"),]$sum_B#9/9 negative
modsums[which(modsums$bestmod=="growm"),]$grow_B#2/2 negative

#Plot a subset of species growth ~dormant season
alldat$sumt<-as.numeric(alldat$sumt)
cedlib<-alldat[alldat$species=="Cedrus_libani",]
pinech<-alldat[alldat$species=="Pinus_echinata",]#
pinfle<-alldat[alldat$species=="Pinus_flexilis",]
chathy<-alldat[alldat$species=="Chamaecyparis_thyoides",]#
pinthu<-alldat[alldat$species=="Pinus_thunbergii",]
thupl<-alldat[alldat$species=="Thuja_plicata",]
picobo<-alldat[alldat$species=="Picea_obovata",]
picpur<-alldat[alldat$species=="Picea_purpurea",]
pinres<-alldat[alldat$species=="Pinus_resinosa",]

colors<-brewer.pal(12,"Paired")

quartz()
plot(c(.5,1),type="p",xlim=c(19,23),ylim=c(.5,2),xlab="Summer mean temperature", ylab="Growth (rwi)",bty="l")
points(thupl$rwi~thupl$sumt, cex=.8, pch=21, bg="lightgray", col="lightgray")
points(chathy$rwi~chathy$sumt, cex=.8, pch=21, bg=colors[7], col=colors[7])
#points(picpur$rwi~picpur$sumt, cex=.8, pch=21, bg=colors[9], col=colors[9])
points(pinech$rwi~pinech$sumt,cex=.8, pch=21,bg=colors[1],col=colors[1])
#points(picobo$rwi~picobo$sumt,cex=.8, pch=21,bg=colors[3],col=colors[3])
points(pinres$rwi~pinres$sumt,cex=.8, pch=21,bg=colors[3],col=colors[3])

#plot lines
sum<-lmer(rwi~sumt+(1|treeid),REML=FALSE,data=pinres)
abline(a=fixef(sum)[1],b=fixef(sum)[2], lwd=2, col=colors[4])
sum<-lmer(rwi~sumt+(1|treeid),REML=FALSE, data=thupl)
abline(a=fixef(sum)[1],b=fixef(sum)[2], lwd=2, col="gray2")
sum<-lmer(rwi~sumt+(1|treeid),REML=FALSE, data=pinech)
abline(a=fixef(sum)[1],b=fixef(sum)[2], lwd=2,  col=colors[2])
#sum<-lm(rwi~sumt,data=picobo)
#abline(a=coef(sum)[1],b=coef(sum)[2], lwd=2, col=colors[4])
sum<-lmer(rwi~sumt+(1|treeid),REML=FALSE,data=chathy)
abline(a=fixef(sum)[1],b=fixef(sum)[2], lwd=2, col=colors[8])
#sum<-lmer(rwi~sumt+(1|treeid),REML=FALSE,data=picpur)
#abline(a=fixef(sum)[1],b=fixef(sum)[2], lwd=2, col=colors[10])

mtext("Chamaecyparis thyoides", col=colors[8])
mtext("Pinus echinata", col=colors[2])
mtext("Thuja plicata", col="gray2")
mtext("Pinus resinosa", col=colors[4])

#sum<-lmer(rwi~sumt+(1|treeid),REML=FALSE,data=picpur)
#abline(a=fixef(sum)[1],b=fixef(sum)[2], lwd=2, col=colors[2])
