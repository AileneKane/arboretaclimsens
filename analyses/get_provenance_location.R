#############################################
### Script to get list of all cored trees ###
### and their population source location  ###
##### (in the form of lat/long or name) #####
######### started April 21, 2026 ############
########### by Ailene Ettinger ##############
######### ailene.ettinger@tnc.org ###########
#############################################

# housekeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# load libraries
library(dplyr)

setwd("~/GitHub/arboretaclimsens")

# read in csv file (exported from Arnold Arboretum's BGBASE in 2012)
bgb <-read.csv("data/ArnArbACCESSIONSALLformetaLAB23Oct2012.csv", header=TRUE)

# read in 2 csv files with ring widths of cored trees 
ang <-read.csv("data/ArnArboringwidths_processed_angiosperms_2018Aug2.csv", header=TRUE)
gym <-read.csv("data/ArnArboringwidths_processed_gymnosperms.csv", header=TRUE)

# combine gym and ang to get list of accession numbers of cored trees

# ugh- cleaning is required! (move the below to a cleaning file sometime)

ang$Site[ang$Site=="Phellodendron_amurense BBM"]<- "Arnold Arboretum-BBM"
ang$Site[ang$Site== "Arn Arboretum" ]<- "Arnold Arboretum"
ang$Site[ang$Site== "Arnold Arb"]<- "Arnold Arboretum"
ang$Site[ang$Site== "populus deltoides"]<- "Arnold Arboretum"
ang$Site[ang$Site== "Phellodendron amurenseBBM"]<- "Arnold Arboretum-BBM"
ang$Site[ang$Site== "bbm"]<- "Arnold Arboretum-BBM"
ang$Site[ang$Site== "Kalopanax septemlobus BBM"]<- "Arnold Arboretum-BBM"
ang$Site[ang$Site== "Kalopanax septemlobus BBM"]<- "Arnold Arboretum-BBM"
ang$Site[ang$Site== "BBM"]<- "Arnold Arboretum-BBM"
ang$Site[ang$Site== "Arb Arboretum" ]<- "Arnold Arboretum"
ang$Site[ang$Site== "Fagus grandifolia" ]<- "Arnold Arboretum"

gymsub<-subset(gym,select=c("TreeName", "Site","Species"))
angsub<-subset(ang,select=c("TreeName", "Site","Species"))

coredtrees<-rbind(gymsub,angsub)
# Select out just the Accession number from the tree name, Accession number is the provenance/source population of the tree
coredtrees$ACC_NUM<- sub("\\*.*$", "", coredtrees$TreeName)

# now pull out location data from bgbase
bgbsub<-subset(bgb, select=c("ACC_NUM","ABBREV_NAME","COUNTRY_FULL","SUB_CNT1","SUB_CNT2","SUB_CNT3","LOCALITY","LAT_DEGREE" ,"LAT_MINUTE","LAT_SECOND" ,"LAT_DIR" ,"LONG_DEGREE","LONG_MINUTE","LONG_SECOND","LONG_DIR","ALTITUDE","ALTITUDE_UNIT","DESCRIPTION","COLLECTION_MISC"))       

# merge the location data into the cored trees data frame
coredtrees_wloc<-left_join(coredtrees,bgbsub, copy=TRUE)

# how many rows are blank for latitude (out of 394 rows total)?

length(which(coredtrees_wloc$LAT_DEGREE==""))
#283

 write.csv(coredtrees_wloc,"output/coredtrees_provenance.csv", row.names=FALSE)      
 