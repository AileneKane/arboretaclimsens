#Cleaning Arnold Arboretum trait data and getting all traits in one file for analyses
#by Ailene
#Started 19 March 2018
## to start
la<-la[1:1819,1:6]
lm<-lm[1:1720,]
#check which ones are missing leaf mass:
unique(lm$AccNum[which(!lm$AccNum %in% la$AccessionNum)])#17 individuals missing la
unique(lm$Species[which(!lm$AccNum %in% la$AccessionNum)])#1 species missing lm
#"CARGLA","CAROVA","KALSEP","PINDEN","PINNIG","PLAOCC","POPDELsspDEL","SORYUA"       "STYJAP"       "CALDEC" 
#which ones are missing leaf area?
unique(la$AccessionNum[which(!la$AccessionNum %in% lm$AccNum)])#17 individuals missing lm
#Fix some typos/data entry mistakes
la$AccessionNum[la$AccessionNum=="224_49_CalDec_2015_001.jpg"]<-"224-49*A"
la$AccessionNum[la$AccessionNum=="1661*K"]<-"16611*K"
la$AccessionNum[la$AccessionNum=="1239-93*C"]<-"1239-83*C"
la$AccessionNum[la$AccessionNum=="17250*B"]<-"1726-81*B"#I'm guessing a little bt here...
la$AccessionNum[la$AccessionNum=="17250*A"]<-"1726-81*A"#I'm guessing a little bt here...
la$AccessionNum[la$AccessionNum=="611-2010*A"]<-"20101*A"#I'm guessing- need to check envelopes##analyzed 7/15/2015 by Silvia, no trees collected with this accession number- must be something else...
la$AccessionNum[la$AccessionNum=="135-38*B"]<-"2019*B"#I'm guessing- need to check envelopes#analyzed 7/15/2015 by Silvia
la$AccessionNum[la$AccessionNum=="135-38*C"]<-"2019*C"#I'm guessing- need to check envelopes#analyzed 7/15/2015 by Silvia
la$AccessionNum[la$AccessionNum=="22886*D"]<-"2019*O"#I'm guessing- need to check envelopes#analyzed 7/15/2015 by Silvia
la$AccessionNum[la$AccessionNum=="21815*E"]<-"2019*I"#I'm guessing- need to check envelopes#analyzed 7/15/2015 by Silvia, no trees collected with this accession number- must be something else...

#The below need to be figured out- not sure what individuals these really are
#la$AccessionNum[la$AccessionNum=="49-2005*A"]#no trees collected with this accession number- must be something else...
#one of the LArix sp? 
#la$AccessionNum[la$AccessionNum=="841-81*A"]#analyzed 7/15/2015 by Silvia and 11/11/15 by sally
#la$AccessionNum[la$AccessionNum=="227-2011*A"]#11/11/15 by sally


lm$AccNum[lm$AccNum=="16526*C"]<-"16536*C"
lm$AccNum[lm$AccNum=="12353*C"]<-"12453*C"
lm$Species[lm$AccNum=="7873*A"]<-"CARGLA"
lm$AccNum[lm$AccNum=="304-87*A"]<-"304-78*A"
lm$AccNum[lm$AccNum=="1539-80*D"]<-"1894-80*D" 
lm$AccNum[lm$AccNum=="1894-80*U"]<-"1894-80*K" 

unique(lm$AccNum[which(!lm$AccNum %in% la$AccessionNum)])#2individuals missing la, for which we have leaf mass
#"20095*A" "7873*A"#these are all from 6/24/2015
unique(la$AccessionNum[which(!la$AccessionNum %in% lm$AccNum)])#8 individuals missing lm
#[1] "49-2005*A"  "611-2010*A" "21815*E"    "22886*D"    "135-38*C"   "135-38*B"  
#[7] "17250*A"    "841-81*A"   "227-2011*A" "17250*B"
colnames(la)[2]<-"AccNum"
lmla<-join(la,lm,by= c("AccNum","LeafIDNum"))
lmla2<-subset(lmla,select=c("AccNum", "LeafIDNum", "leaf_area_cm","leaf_mass_mg"))
lmla2$lma<-as.numeric(lmla2$leaf_mass_mg)/as.numeric(lmla2$leaf_area_cm)
lma.mn<-tapply(lmla2$lma,lmla2$AccNum, mean, na.rm=T)
lma<-as.data.frame(cbind(names(lma.mn),lma.mn))
colnames(lma)[1]<-"AccNum"
rownames(lma)<-NULL
lma$lma.mn<-as.numeric(lma$lma.mn)
lma<-na.omit(lma)
woodmass<-woodmass[,1:11]
woodmass$dens<-(woodmass$TwigDryMass_mg/1000)/woodmass$TwigVolume_cm3
wdens<-subset(woodmass,select=c("Species","AccessionNumber", "dens"))
colnames(wdens)[2]<-"AccNum"
wdens$dens[wdens$dens<0]<-abs(wdens$dens[wdens$dens<0])#remove negative dens and replace with positive
wdens$dens[wdens$dens>1]<-wdens$dens[wdens$dens>1]-1
lma.wdens<-join(lma,wdens,by="AccNum")#303 individuals for which we have both lma and wood density
#now add height
height$height1<-as.numeric(height$height1)
height$height2<-as.numeric(height$height2)
height$mn.height<-rowMeans(height[,8:9],na.rm = FALSE)#could look at arboretum data to get missing values
height2<-subset(height,select=c("accession_no","genus","species","var_ssp","mn.height"))
colnames(height2)[1]<-"AccNum"
lma.wdens.ht<-merge(lma.wdens,height2,by="AccNum")#319 individuals for which we have all trait and height data
head(lma.wdens.ht)
lma.wdens.ht[lma.wdens.ht=="<NA>"] =NA
lma.wdens.ht[lma.wdens.ht=="NaN"] =NA
lma.wdens.ht$lma.mn<-as.numeric(as.character(lma.wdens.ht$lma.mn))
lma.wdens.ht$dens<-as.numeric(as.character(lma.wdens.ht$dens))
lma.wdens.ht$mn.height<-as.numeric(as.character(lma.wdens.ht$mn.height))
lma.wdens.ht$SpName<-paste(lma.wdens.ht$genus,lma.wdens.ht$species,sep = "_")
#add angio/gymnosperm
lma.wdens.ht$group<-NA
lma.wdens.ht[lma.wdens.ht$genus=="Tsuga"|lma.wdens.ht$genus=="Pinus"|lma.wdens.ht$genus=="Larix"|lma.wdens.ht$genus=="Cryptomeria"|lma.wdens.ht$genus=="Thuja"|lma.wdens.ht$genus=="Pseudotsuga"|lma.wdens.ht$genus=="Chamaecyparis"|lma.wdens.ht$genus=="Picea"|lma.wdens.ht$genus=="Abies"|lma.wdens.ht$genus=="Taxus"|lma.wdens.ht$genus=="Cedrus"|lma.wdens.ht$genus=="Calocedrus"|lma.wdens.ht$genus=="Metasequoia",]$group<-"gymnosperm"
lma.wdens.ht[lma.wdens.ht$genus=="Liriodendron"|lma.wdens.ht$genus=="Pyrus"|lma.wdens.ht$genus=="Phellodendron"|lma.wdens.ht$genus=="Platanus"|lma.wdens.ht$genus=="Crataegus"|lma.wdens.ht$genus=="Betula"|lma.wdens.ht$genus=="Quercus"|lma.wdens.ht$genus=="Styphnolobium"|lma.wdens.ht$genus=="Gleditsia"|lma.wdens.ht$genus=="Catalpa"|lma.wdens.ht$genus=="Kalopanax"|lma.wdens.ht$genus=="Liquidambar"|lma.wdens.ht$genus=="Aesculus"|lma.wdens.ht$genus=="Carya"|lma.wdens.ht$genus=="Tilia"|lma.wdens.ht$genus=="Fagus"|lma.wdens.ht$genus=="Fraxinus"|lma.wdens.ht$genus=="Sorbus"|lma.wdens.ht$genus=="Populus"|lma.wdens.ht$genus=="Ulmus"|lma.wdens.ht$genus=="Zelkova",]$group<-"angiosperm"
#Merge spring phenology with fall phenology; and (separately) with other traits, by individual
springphen2<-springphen[,3:8]
colnames(springphen2)[1:2]<-c("AccNum","SpName")
fallphen<-fallphen[,3:13]
colnames(fallphen)[1:2]<-c("AccNum","SpName")
fallphen2<-subset(fallphen,select=c("AccNum","FirstColDOY", "FirstDropDOY", "PeakDropDOY"))
allphen<-inner_join(fallphen2,springphen2,by="AccNum")
colnames(allphen)[1:5]<-c("accnum","col","drop","pkdrop","sp")
colnames(allphen)[8:9]<-c("bb","le")

traits<-left_join(springphen2,lma.wdens.ht)
traits$SpName[traits$SpName=="Pinusflexilis"]<-"Pinus_flexilis"
traits$SpName[traits$SpName=="Carpinuscaroliniana"]<-"Carpinus_caroliniana"
traits<-subset(traits,select=c("AccNum", "SpName", "genus","group","LeafBBdoy","LeafEmergdoy","lma.mn","dens","mn.height"))
colnames(traits)[1:2]<-c("accnum","sp")
colnames(traits)[5]<-"bb"
colnames(traits)[6]<-"le"
colnames(traits)[7]<-"lma"
colnames(traits)[9]<-"ht"

