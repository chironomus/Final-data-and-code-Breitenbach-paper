#starting from the original excel files, only modification is to remove species with  low dance manually, #since it is faster in excel, all other manipulations 
#done on this original data via code
#removed spp had cumulative abundance (41 yrs sum <1000 specimens total)

# subsetting and preparing Breitenbach datasets for the analysis#
#package  installation
#install.packages(c("plyr", "dplyr", "reshape2", "ggplot2", "nlme","gridExtra",
 "FD","mblm","trend","broom","ggfortify","graphics","Kendall","DataCombine",
 "modifiedmk","boot","vegan","gtools","codyn"))
require(plyr)
require(dplyr)
require(forecast)
require(tidyverse)
require(reshape2)
require(nlme)
require(ggplot2)
require(gridExtra)
require(FD)
require(mblm)
require(trend)
require(broom)
require(ggfortify)
library(DataCombine)
require(graphics)
require(Kendall)
require(boot)
require(modifiedmk)
require(vegan)
library(gtools)
library(vegan)
###########################################
##########################################
#pnenological and abundance trends of the 31 most abundant taxa analyzed

#set your working directory with uploaded data
#setwd()
#this is just my directory, plug yours here
setwd("C:/Users/Viktor/Documents/Baranov/my papers/Submitted/4_Breitenbach Emergenzdaten/portable data and code")
 
############################################################
############################################################
#read the data for the Trap B post 2005
bp2005= read.csv("siteBpost20051.csv",sep=";")

head(bp2005)

#replacing NA by "0" since NA in this matrix is a zero abundance of the taxa

bp2005[is.na(bp2005)] <- 0 

#subset taxa
taxa=bp2005[,4:65]
# summarizing males and females into single sp column
output <-data.frame( Bae_ver=apply(taxa[1:2], 1, sum),Bae_roh=apply(taxa[3:4], 1, sum), Eph_muc=apply(taxa[5:6], 1, sum),Eph_ign=apply(taxa[7:8],1, sum),
                     Cen_lut=apply(taxa[9:10],1, sum),Par_sub=apply(taxa[11:12],1, sum),Amp_sta=apply(taxa[13:14],1, sum),Bra_ris=apply(taxa[15:16],1, sum),
                     Iso_goe=apply(taxa[17:18],1, sum),Leu_dig=apply(taxa[19:20],1, sum),Leu_nig=apply(taxa[21:22],1, sum), Leu_pri=apply(taxa[23:24],1, sum),
                     Nem_cam=apply(taxa[25:26],1, sum),Nem_fle=apply(taxa[27:28],1, sum),Nem_mar=apply(taxa[29:30],1, sum),Nrl_pic=apply(taxa[31:32],1, sum),
                     Pro_aub=apply(taxa[33:34],1, sum),Pro_int=apply(taxa[35:36],1, sum),Pro_mey=apply(taxa[37:38],1, sum),Sip_tor=apply(taxa[39:40],1, sum),
                     Aga_fus=apply(taxa[41:42],1, sum),Apa_fim=apply(taxa[43:44],1, sum),Cha_vil=apply(taxa[45:46],1, sum),Dru_ann=apply(taxa[47:48],1, sum),
                     Ple_con=apply(taxa[49:50],1, sum),Pot_luc=apply(taxa[51:52],1, sum),Rhy_fas=apply(taxa[53:54],1, sum),Ser_per=apply(taxa[55:56],1, sum),
                     Sil_pal=apply(taxa[57:58],1, sum),Tin_ros=apply(taxa[59:60],1, sum),Wor_occ=apply(taxa[61:62],1, sum))

bp2005=cbind(bp2005,output)# combine column sum of male+female abundance with date data
bp2005=as.data.frame(bp2005)#as data frame
#adding datein r readible date format
bp2005$date <- as.Date(with(bp2005, paste(B_Jahr,B_Monat,B_Tag,sep="-")), "%Y-%m-%d")
#calculating sequential number of week (in a year) for each collection date
weeknum=strftime(bp2005$date, format = "%V")
bp2005$weeknum=as.numeric(weeknum)#attach week number id to the dataset as a column
# add dummy "iday" variable for compatibility with the other datasets, "iday" is a number of the days in the year, used in some datasets

iday=bp2005$weeknum*7

bp2005$iday=iday

#reorder columns to match other datasets
bp2005=bp2005[c(1:3,99,97,98,4:96)]
############################################################
############################################################
#Same procedures as in lines 37-71 are repeated for the dataset from the trap B for years 1986-2006
#preparation of the dataset for site B 1982-2006
b82_2006= read.csv("siteBpost82_2006.csv",sep=";")

#replacing NA by "0" since NA in this matrix is a zero abundance of the taxa
b82_2006[is.na(b82_2006)] <- 0
#subset taxa
taxa=b82_2006[,7:66]
# summarizing males and females into single sp column ! NOTE Bae_ver variable is invariably identified as factor by r, so I have manually calculated m+f sum for this 
#and add it to the table
output <-data.frame(Bae_roh=apply(taxa[1:2], 1, sum), Eph_muc=apply(taxa[3:4], 1, sum),Eph_ign=apply(taxa[5:6],1, sum),
                    Cen_lut=apply(taxa[7:8],1, sum),Par_sub=apply(taxa[9:10],1, sum),Amp_sta=apply(taxa[11:12],1, sum),Bra_ris=apply(taxa[13:14],1, sum),
                    Iso_goe=apply(taxa[15:16],1, sum),Leu_dig=apply(taxa[17:18],1, sum),Leu_nig=apply(taxa[19:20],1, sum), Leu_pri=apply(taxa[21:22],1, sum),
                    Nem_cam=apply(taxa[23:24],1, sum),Nem_fle=apply(taxa[25:26],1, sum),Nem_mar=apply(taxa[27:28],1, sum),Nrl_pic=apply(taxa[29:30],1, sum),
                    Pro_aub=apply(taxa[31:32],1, sum),Pro_int=apply(taxa[33:34],1, sum),Pro_mey=apply(taxa[35:36],1, sum),Sip_tor=apply(taxa[37:38],1, sum),
                    Aga_fus=apply(taxa[39:40],1, sum),Apa_fim=apply(taxa[41:42],1, sum),Cha_vil=apply(taxa[43:44],1, sum),Dru_ann=apply(taxa[45:46],1, sum),
                    Ple_con=apply(taxa[47:48],1, sum),Pot_luc=apply(taxa[49:50],1, sum),Rhy_fas=apply(taxa[51:52],1, sum),Ser_per=apply(taxa[53:54],1, sum),
                    Sil_pal=apply(taxa[55:56],1, sum),Tin_ros=apply(taxa[57:58],1, sum),Wor_occ=apply(taxa[59:60],1, sum))

b82_2006=cbind(b82_2006,output)
b82_2006=as.data.frame(b82_2006)
#adding date
b82_2006$date <- as.Date(with(b82_2006, paste(B_Jahr,B_Monat,B_Tag,sep="-")), "%Y-%m-%d")
#calculating number of week (in a year) for each collection date
weeknum=strftime(b82_2006$date, format = "%V")
b82_2006$weeknum=as.numeric(weeknum)

#reorder columns to match other datasets
b82_2006=b82_2006[c(1:4,98,99,5:97)]


############################################################
############################################################
#same procedures repeared for the trap B 1971-1972 preparation

b71= read.csv("siteB71_82.csv",sep=";")
#replacing NA by "0" since NA in this matrix is a zero abundance of the taxa
b71[is.na(b71)] <- 0
#subset taxa
taxa=b71[,5:66] 
# summarizing males and females into single sp column
output <-data.frame( Bae_ver=apply(taxa[1:2], 1, sum),Bae_roh=apply(taxa[3:4], 1, sum), Eph_muc=apply(taxa[5:6], 1, sum),Eph_ign=apply(taxa[7:8],1, sum),
                     Cen_lut=apply(taxa[9:10],1, sum),Par_sub=apply(taxa[11:12],1, sum),Amp_sta=apply(taxa[13:14],1, sum),Bra_ris=apply(taxa[15:16],1, sum),
                     Iso_goe=apply(taxa[17:18],1, sum),Leu_dig=apply(taxa[19:20],1, sum),Leu_nig=apply(taxa[21:22],1, sum), Leu_pri=apply(taxa[23:24],1, sum),
                     Nem_cam=apply(taxa[25:26],1, sum),Nem_fle=apply(taxa[27:28],1, sum),Nem_mar=apply(taxa[29:30],1, sum),Nrl_pic=apply(taxa[31:32],1, sum),
                     Pro_aub=apply(taxa[33:34],1, sum),Pro_int=apply(taxa[35:36],1, sum),Pro_mey=apply(taxa[37:38],1, sum),Sip_tor=apply(taxa[39:40],1, sum),
                     Aga_fus=apply(taxa[41:42],1, sum),Apa_fim=apply(taxa[43:44],1, sum),Cha_vil=apply(taxa[45:46],1, sum),Dru_ann=apply(taxa[47:48],1, sum),
                     Ple_con=apply(taxa[49:50],1, sum),Pot_luc=apply(taxa[51:52],1, sum),Rhy_fas=apply(taxa[53:54],1, sum),Ser_per=apply(taxa[55:56],1, sum),
                     Sil_pal=apply(taxa[57:58],1, sum),Tin_ros=apply(taxa[59:60],1, sum),Wor_occ=apply(taxa[61:62],1, sum))

b71=cbind(b71,output)
b71=as.data.frame(b71)
#adding date
b71$date <- as.Date(with(b71, paste(B_Jahr,B_Monat,B_Tag,sep="-")), "%Y-%m-%d")
#calculating number of week (in a year) for each collection date
weeknum=strftime(b71$date, format = "%V")
b71$weeknum=as.numeric(weeknum)

#reorder columns to match other datasets
b71=b71[c(1:4,98,99,5:97)]
############################################################
############################################################
#same procedures repeared for the trap B 1969-1971 preparation

# site B1969-1971 preparaiton
b69= read.csv("siteB69_71.csv",sep=";")
#replacing NA by "0" since NA in this matrix is a zero abundance of the taxa
b69[is.na(b69)] <- 0
#subset taxa
taxa=b69[,5:66] 
# summarizing males and females into single sp column
output <-data.frame( Bae_ver=apply(taxa[1:2], 1, sum),Bae_roh=apply(taxa[3:4], 1, sum), Eph_muc=apply(taxa[5:6], 1, sum),Eph_ign=apply(taxa[7:8],1, sum),
                     Cen_lut=apply(taxa[9:10],1, sum),Par_sub=apply(taxa[11:12],1, sum),Amp_sta=apply(taxa[13:14],1, sum),Bra_ris=apply(taxa[15:16],1, sum),
                     Iso_goe=apply(taxa[17:18],1, sum),Leu_dig=apply(taxa[19:20],1, sum),Leu_nig=apply(taxa[21:22],1, sum), Leu_pri=apply(taxa[23:24],1, sum),
                     Nem_cam=apply(taxa[25:26],1, sum),Nem_fle=apply(taxa[27:28],1, sum),Nem_mar=apply(taxa[29:30],1, sum),Nrl_pic=apply(taxa[31:32],1, sum),
                     Pro_aub=apply(taxa[33:34],1, sum),Pro_int=apply(taxa[35:36],1, sum),Pro_mey=apply(taxa[37:38],1, sum),Sip_tor=apply(taxa[39:40],1, sum),
                     Aga_fus=apply(taxa[41:42],1, sum),Apa_fim=apply(taxa[43:44],1, sum),Cha_vil=apply(taxa[45:46],1, sum),Dru_ann=apply(taxa[47:48],1, sum),
                     Ple_con=apply(taxa[49:50],1, sum),Pot_luc=apply(taxa[51:52],1, sum),Rhy_fas=apply(taxa[53:54],1, sum),Ser_per=apply(taxa[55:56],1, sum),
                     Sil_pal=apply(taxa[57:58],1, sum),Tin_ros=apply(taxa[59:60],1, sum),Wor_occ=apply(taxa[61:62],1, sum))
b69=cbind(b69,output)
b69=as.data.frame(b69)
#adding date
b69$date <- as.Date(with(b69, paste(B_Jahr,B_Monat,B_Tag,sep="-")), "%Y-%m-%d")
#calculating number of week (in a year) for each collection date
weeknum=strftime(b69$date, format = "%V")
b69$weeknum=as.numeric(weeknum)

#reorder columns to match other datasets
b69=b69[c(1:4,98,99,5:97)]


############################################################
############################################################
#compile single dataframe for all the dataframes covering the trap B (1969-2010)

b1=rbind(b69,b71)
b2=rbind(b1,b82_2006)
siteB_whole=rbind(b2,bp2005) # mind 2006 overlap between site B 82-06 and "post 2005" data NOTe. data are identical, 2006 deleted from "post 2005" datset
#Wagner 3.11.2017, pers.comm

#subset species level sums for males and females, delete males and females individual columns for now

ept2<- siteB_whole[c(1:6,69:99)]

ept2<- subset(ept2, B_Jahr != "0") #drop abberant  "0" from the "year" column

ept2=ept2[order(ept2$B_Jahr,ept2$iday),]#order by the year, number of day
#split dataframe into list by year
mylist <- split(ept2,ept2$B_Jahr) 

#each year is now separated as an object within list "mylist"
#that allows us to calculate phenological variables+abundance for each spp for each year

#1969
e69=as.data.frame(mylist[[1]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e69, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))

xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]
# central tendency of phenology according to  Ash 2015 doi:  10.1073/pnas.1421946112
cumper1969<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))

#date of first emergence per sp per year 
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first69<-as.data.frame(do.call(rbind,first_em))

#date of last emergence per sp per year 
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last69<-as.data.frame(do.call(rbind,last_em))

###########################################
#Identical procedure for all the proceeding years
#1970
e70=as.data.frame(mylist[[2]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e70, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1970<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first70<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last70<-as.data.frame(do.call(rbind,last_em))


#1971
e71=as.data.frame(mylist[[3]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e71, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1971<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first71<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last71<-as.data.frame(do.call(rbind,last_em))


#1972
e72=as.data.frame(mylist[[4]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e72, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1972<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first72<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last72<-as.data.frame(do.call(rbind,last_em))

#1973
e73=as.data.frame(mylist[[5]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e73, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1973<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first73<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last73<-as.data.frame(do.call(rbind,last_em))

#1974
e74=as.data.frame(mylist[[6]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e74, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1974<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first74<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last74<-as.data.frame(do.call(rbind,last_em))


#1975
e75=as.data.frame(mylist[[7]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e75, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1975<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first75<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last75<-as.data.frame(do.call(rbind,last_em))


#1976
e76=as.data.frame(mylist[[8]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e76, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1976<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first76<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last76<-as.data.frame(do.call(rbind,last_em))



#1977
e77=as.data.frame(mylist[[9]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e77, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1977<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015
#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first77<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last77<-as.data.frame(do.call(rbind,last_em))

#1978
e78=as.data.frame(mylist[[10]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e78, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1978<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first78<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last78<-as.data.frame(do.call(rbind,last_em))



#1979
e79=as.data.frame(mylist[[11]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e79, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1979<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first79<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last79<-as.data.frame(do.call(rbind,last_em))



#1980
e80=as.data.frame(mylist[[12]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e80, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1980<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first80<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last80<-as.data.frame(do.call(rbind,last_em))



#1981
e81=as.data.frame(mylist[[13]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e81, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1981<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first81<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last81<-as.data.frame(do.call(rbind,last_em))

#1982
e82=as.data.frame(mylist[[14]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e82, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1982<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first82<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last82<-as.data.frame(do.call(rbind,last_em))

#1983
e83=as.data.frame(mylist[[15]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e83, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1983<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first83<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last83<-as.data.frame(do.call(rbind,last_em))


#1984
e84=as.data.frame(mylist[[16]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e84, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1984<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first84<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last84<-as.data.frame(do.call(rbind,last_em))

#1985
e85=as.data.frame(mylist[[17]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e85, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1985<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first85<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last85<-as.data.frame(do.call(rbind,last_em))



#1986
e86=as.data.frame(mylist[[18]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e86, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1986<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first86<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last86<-as.data.frame(do.call(rbind,last_em))



#1987
e87=as.data.frame(mylist[[19]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e87, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1987<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first87<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last87<-as.data.frame(do.call(rbind,last_em))


#1988
e88=as.data.frame(mylist[[20]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e88, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1988<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first88<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last88<-as.data.frame(do.call(rbind,last_em))


#1989
e89=as.data.frame(mylist[[21]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e89, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1989<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first89<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last89<-as.data.frame(do.call(rbind,last_em))



#1990
e90=as.data.frame(mylist[[22]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e90, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1990<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first90<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last90<-as.data.frame(do.call(rbind,last_em))



#1991
e91=as.data.frame(mylist[[23]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e91, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1991<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first91<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last91<-as.data.frame(do.call(rbind,last_em))



#1992
e92=as.data.frame(mylist[[24]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e92, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1992<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first92<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last92<-as.data.frame(do.call(rbind,last_em))



#1993
e93=as.data.frame(mylist[[25]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e93, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1993<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first93<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last93<-as.data.frame(do.call(rbind,last_em))


#1994
e94=as.data.frame(mylist[[26]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e94, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1994<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first94<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last94<-as.data.frame(do.call(rbind,last_em))



#1995
e95=as.data.frame(mylist[[27]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e95, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1995<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first95<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last95<-as.data.frame(do.call(rbind,last_em))



#1996
e96=as.data.frame(mylist[[28]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e96, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1996<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first96<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last96<-as.data.frame(do.call(rbind,last_em))


#1997
e97=as.data.frame(mylist[[29]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e97, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1997<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first97<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last97<-as.data.frame(do.call(rbind,last_em))

#1998
e98=as.data.frame(mylist[[30]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e98, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1998<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first98<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last98<-as.data.frame(do.call(rbind,last_em))


#1999
e99=as.data.frame(mylist[[31]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e99, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper1999<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first99<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last99<-as.data.frame(do.call(rbind,last_em))


#2000
e00=as.data.frame(mylist[[32]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e00, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper2000<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first00<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last00<-as.data.frame(do.call(rbind,last_em))



#2001
e01=as.data.frame(mylist[[33]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e01, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper2001<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first01<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last01<-as.data.frame(do.call(rbind,last_em))


#2002
e02=as.data.frame(mylist[[34]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e02, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper2002<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2025


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first02<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last02<-as.data.frame(do.call(rbind,last_em))



#2003
e03=as.data.frame(mylist[[35]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e03, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper2003<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2035


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first03<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last03<-as.data.frame(do.call(rbind,last_em))



#2004
e04=as.data.frame(mylist[[36]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e04, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper2004<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2045


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first04<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last04<-as.data.frame(do.call(rbind,last_em))


#2005
e05=as.data.frame(mylist[[37]])
#melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e05, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper2005<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2055


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first05<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last05<-as.data.frame(do.call(rbind,last_em))



#2006
e06=as.data.frame(mylist[[38]])
melt dataframe to long matrix with all taxa in column "variable" and all abundances in column "value"
eptB_melt <- melt(e06, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
xy<-eptB_melt[!(apply(eptB_melt,1,function(y) any(y==0))),]

cumper2006<-ddply(xy,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2065


#first emergence day
first_em<-lapply(unique(xy$variable),function(a)head (subset(xy,variable==a & value>0),1))
first06<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(xy$variable),function(a)tail (subset(xy,variable==a & value>0),1))
last06<-as.data.frame(do.call(rbind,last_em))

# pre and post 2006 datasets should first be combined seprately, since 1969-2006 requiring division of all phenological variables values  by 7 to get week numbers

#2007
e07=as.data.frame(mylist[[39]])
eptB_melt <-  melt(e07, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))


cumper2007<-ddply(eptB_melt,.(variable),transform,ct=sum(value*weeknum)/sum(value))# central tendency according to  Ash 2015

#first emergence day
first_em<-lapply(unique(eptB_melt$variable),function(a)head (subset(eptB_melt,variable==a & value>0),1))
first07<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(eptB_melt$variable),function(a)tail (subset(eptB_melt,variable==a & value>0),1))
last07<-as.data.frame(do.call(rbind,last_em))

#2008
e08=as.data.frame(mylist[[40]])
eptB_melt <-  melt(e08, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
eptB_melt=na.omit(eptB_melt)
cumper2008<-ddply(eptB_melt,.(variable),transform,ct=sum(value*weeknum)/sum(value))# central tendency according to  Ash 2015

#first emergence dayu
first_em<-lapply(unique(eptB_melt$variable),function(a)head (subset(eptB_melt,variable==a & value>0),1))
first08<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(eptB_melt$variable),function(a)tail (subset(eptB_melt,variable==a & value>0),1))
last08<-as.data.frame(do.call(rbind,last_em))
#2009
e09=as.data.frame(mylist[[41]])
eptB_melt <-  melt(e09, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))


cumper2009<-ddply(eptB_melt,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015

#first emergence day
first_em<-lapply(unique(eptB_melt$variable),function(a)head (subset(eptB_melt,variable==a & value>0),1))
first09<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(eptB_melt$variable),function(a)tail (subset(eptB_melt,variable==a & value>0),1))
last09<-as.data.frame(do.call(rbind,last_em))


#2010
e10=as.data.frame(mylist[[42]])
eptB_melt <-  melt(e10, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))


cumper2010<-ddply(eptB_melt,.(variable),transform,ct=sum(value*iday)/sum(value))# central tendency according to  Ash 2015

#first emergence day
first_em<-lapply(unique(eptB_melt$variable),function(a)head (subset(eptB_melt,variable==a & value>0),1))
first10<-as.data.frame(do.call(rbind,first_em))

#last emergence day
last_em<-lapply(unique(eptB_melt$variable),function(a)tail (subset(eptB_melt,variable==a & value>0),1))
last10<-as.data.frame(do.call(rbind,last_em))

#merge phenological variables for each year into a dataframes, 
#which in turn would be merged in the final dataframe
# I switched off writing of the each datframe into working directory, 
#since its slowing computation considerably
ept1969<-merge(cumper1969, merge(first69,last69, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1969, file = "eptphen1969.csv")

ept1970<-merge(cumper1970, merge(first70,last70, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1970, file = "eptphen1970.csv")

ept1971<-merge(cumper1971, merge(first71,last71, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1971, file = "eptphen1971.csv")

ept1972<-merge(cumper1972, merge(first72,last72, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1972, file = "eptphen1972.csv")

ept1973<-merge(cumper1973, merge(first73,last73, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1973, file = "eptphen1973.csv")

ept1974<-merge(cumper1974, merge(first74,last74, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1974, file = "eptphen1974.csv")

ept1975<-merge(cumper1975, merge(first75,last75, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1975, file = "eptphen1975.csv")

ept1976<-merge(cumper1976, merge(first76,last76, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1976, file = "eptphen1976.csv")

ept1977<-merge(cumper1977, merge(first77,last77, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1977, file = "eptphen1977.csv")

ept1978<-merge(cumper1978, merge(first78,last78, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1978, file = "eptphen1978.csv")

ept1979<-merge(cumper1979, merge(first79,last79, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1979, file = "eptphen1979.csv")

ept1980<-merge(cumper1980, merge(first80,last80, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1980, file = "eptphen1980.csv")

ept1981<-merge(cumper1981, merge(first81,last81, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1981, file = "eptphen1981.csv")

ept1982<-merge(cumper1982, merge(first82,last82, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1982, file = "eptphen1982.csv")

ept1983<-merge(cumper1983, merge(first83,last83, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1983, file = "eptphen1983.csv")


ept1984<-merge(cumper1984, merge(first84,last84, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1984, file = "eptphen1984.csv")


ept1985<-merge(cumper1985, merge(first85,last85, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1985, file = "eptphen1985.csv")

ept1986<-merge(cumper1986, merge(first86,last86, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1986, file = "eptphen1986.csv")

ept1987<-merge(cumper1987, merge(first87,last87, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1987, file = "eptphen1987.csv")

ept1988<-merge(cumper1988, merge(first88,last88, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1988, file = "eptphen1988.csv")

ept1989<-merge(cumper1989, merge(first89,last89, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1989, file = "eptphen1989.csv")

ept1990<-merge(cumper1990, merge(first90,last90, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1990, file = "eptphen1990.csv")

ept1991<-merge(cumper1991, merge(first91,last91, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1991, file = "eptphen1991.csv")

ept1992<-merge(cumper1992, merge(first92,last92, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1992, file = "eptphen1992.csv")

ept1993<-merge(cumper1993, merge(first93,last93, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1993, file = "eptphen1993.csv")

ept1994<-merge(cumper1994, merge(first94,last94, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1994, file = "eptphen1994.csv")

ept1995<-merge(cumper1995, merge(first95,last95, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1995, file = "eptphen1995.csv")

ept1996<-merge(cumper1996, merge(first96,last96, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1996, file = "eptphen1996.csv")

ept1997<-merge(cumper1997, merge(first97,last97, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1997, file = "eptphen1997.csv")

ept1998<-merge(cumper1998, merge(first98,last98, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1998, file = "eptphen1998.csv")

ept1999<-merge(cumper1999, merge(first99,last99, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept1999, file = "eptphen1999.csv")

ept2000<-merge(cumper2000, merge(first00,last00, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept2000, file = "eptphen2000.csv")

ept2001<-merge(cumper2001, merge(first01,last01, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept2001, file = "eptphen2001.csv")

ept2002<-merge(cumper2002, merge(first02,last02, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept2002, file = "eptphen2002.csv")

ept2003<-merge(cumper2003, merge(first03,last03, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept2003, file = "eptphen2003.csv")

ept2004<-merge(cumper2004, merge(first04,last04, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept2004, file = "eptphen2004.csv")

ept2005<-merge(cumper2005, merge(first05,last05, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept2005, file = "eptphen2005.csv")

ept2006<-merge(cumper2006, merge(first06,last06, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
#write.table(ept2006, file = "eptphen2006.csv")


ept_1<-rbind(ept1969,ept1970,ept1971,ept1972,ept1973,ept1974,ept1975,ept1976,ept1977,ept1978,ept1979,ept1980,ept1981,ept1982,ept1983,ept1984,ept1985,ept1986,ept1987)
ept_2<-rbind(ept1988,ept1989,ept1990,ept1991,ept1992,ept1993,ept1994,ept1995,ept1996,ept1997,ept1998,ept1999,ept2000,ept2001,ept2002,ept2003,ept2004,ept2005, ept2006)
ept_tot<-rbind(ept_1,ept_2)#1969-2006 pheno dataframe
#will rename later ct is central tendency. iday.x=first emergence, iday.y=last emergence
ept_tot$week.first=(ept_tot$iday.x)/7#multiply days of phenophase by 7 to get week numbers
ept_tot$week.last=(ept_tot$iday.y)/7#multiply days of phenophase by 7 to get week numbers
ept_tot$ct.week=(ept_tot$ct)/7#multiply days of phenophase by 7 to get week numbers
#write.table(ept_tot, "siteB1969_2006_phenovars.txt", sep="\t")


#merge 2007-2010 trap B
ept2007<-merge(cumper2007, merge(first07,last07, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
ept2007<-na.omit(ept2007)
ept2007$week.first=(ept2007$iday.x)
ept2007$week.last=(ept2007$iday.y)
ept2007$ct.week=(ept2007$ct)
ept2007$duration=(ept2007$week.last-ept2007$week.first)/7
ept2008<-merge(cumper2008, merge(first08,last08, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
ept2008<-na.omit(ept2008)
ept2008$week.first=(ept2008$iday.x)
ept2008$week.last=(ept2008$iday.y)
ept2008$ct.week=(ept2008$ct)
ept2008$duration=(ept2008$week.last-ept2008$week.first)/7
ept2009<-merge(cumper2009, merge(first09,last09, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
ept2009$week.first=(ept2009$iday.x/7)
ept2009$week.last=(ept2009$iday.y/7)
ept2009$ct.week=(ept2009$ct/7)
ept2009$duration=(ept2009$week.last-ept2009$week.first)
ept2010<-merge(cumper2010, merge(first10,last10, by="variable", all.x=TRUE, all.y=TRUE), by = "variable", all.x = TRUE, all.y = TRUE)
ept2010$week.first=(ept2010$iday.x/7)
ept2010$week.last=(ept2010$iday.y/7)
ept2010$ct.week=(ept2010$ct/7)
ept2010$duration=(ept2010$week.last-ept2010$week.first)

ept_new<-rbind(ept2007,ept2008,ept2009,ept2010) 
#since here all the data has weekly resolution, I don't have to / them by 7 to get phenoweeks, but I will ddivide them by 1 and add the column, in order to match 
#num of cols in 1969-2006 data
#ept_new$week.first=(ept_new$iday.x/)
#ept_new$week.last=(ept_new$iday.y)
#ept_new$ct.week=(ept_new$ct)
ept_tot$duration=ept_tot$week.last-ept_tot$week.first
b1=rbind(ept_tot,ept_new)#total dataset for the phenological analysis
#calculate duration of emergence for every spp
#for every year from subtractin date of first emergence
#from the date of last emergence


b1$Jahr=b1$B_Jahr#rename some variables for convineince
b1$Monat=b1$B_Monat
b1$Tag=b1$B_Tag
s_pheno=b1

#calculate cumulative abundance of every species per year
eptmelt2<- melt(ept2, id=c("B_Jahr","B_Monat","B_Tag","iday","date","weeknum"))
bx1=aggregate(eptmelt2[,8], list(eptmelt2$B_Jahr,eptmelt2$variable), sum) #aggregate means per sp per year fro all 4 phenophases
bx1$variable=bx1$Group.2
#aggregate phenological variables values per sp per year
d1=aggregate(b1[, c(24:27)], list(b1$Jahr,b1$variable), mean) #aggregate means per sp per year fro all 4 phenophases
#merge phenological data and cumulative abundance
d1<-merge(d1,bx1,by =c("Group.2","Group.1"), all.x = TRUE, all.y = TRUE)
d1$value=d1$x #rename
#remove 2006 since year was undersampled
d1=subset(d1,d1$Group.1!="2006")
d1$taxa=d1$Group.2#rename column renamed previously by aggregation
d1$Year=d1$Group.1#rename column renamed previously by aggregation
d1=na.omit(d1)
d1=d1[,c(3:6,9:11)]


#separate 31 most abundant species into dataset
#see supplement containing abbreviations and full names of all spp
#collected
Aga=subset(d1,d1$taxa=="Aga_fus")

Amp=subset(d1,d1$taxa=="Amp_sta")

Apa=subset(d1,d1$taxa=="Apa_fim")


Bro=subset(d1,d1$taxa=="Bae_roh")

Bve=subset(d1,d1$taxa=="Bae_ver")

Bra=subset(d1,d1$taxa=="Bra_ris")

Cen=subset(d1,d1$taxa=="Cen_lut")

Cha=subset(d1,d1$taxa=="Cha_vil")

Dru=subset(d1,d1$taxa=="Dru_ann")

Epi=subset(d1,d1$taxa=="Eph_ign")

Epm=subset(d1,d1$taxa=="Eph_muc")

Iso=subset(d1,d1$taxa=="Iso_goe")

Leu_dig=subset(d1,d1$taxa=="Leu_dig")

Leu_nig=subset(d1,d1$taxa=="Leu_nig")

Leu_pri=subset(d1,d1$taxa=="Leu_pri")

Nem_cam=subset(d1,d1$taxa=="Nem_cam")



Nem_fle=subset(d1,d1$taxa=="Nem_fle")

Nem_mar=subset(d1,d1$taxa=="Nem_mar")

Nrl_pic=subset(d1,d1$taxa=="Nrl_pic")

Par_sub=subset(d1,d1$taxa=="Par_sub")

Ple_con=subset(d1,d1$taxa=="Ple_con")

Pot_luc=subset(d1,d1$taxa=="Pot_luc")

Pro_aub=subset(d1,d1$taxa=="Pro_aub")

Pro_int=subset(d1,d1$taxa=="Pro_int")

Pro_mey=subset(d1,d1$taxa=="Pro_mey")

Rhy_fas=subset(d1,d1$taxa=="Rhy_fas")

Ser_per=subset(d1,d1$taxa=="Ser_per")

Sil_pal=subset(d1,d1$taxa=="Sil_pal")

Sip_tor=subset(d1,d1$taxa=="Sip_tor")

Tin_ros=subset(d1,d1$taxa=="Tin_ros")

Wor_occ=subset(d1,d1$taxa=="Wor_occ")

#run Mann Kendall-tests and Sen's slope models for each of 31 spp
#central tendency trends
mk_aga_ct=tidy(mk.test(Aga$ct.week))
tidy_aga_ct=tidy(mblm(ct.week~Year,data=Aga))
#date of first emergence trends
mk_aga_week.first=tidy(mk.test(Aga$week.first))
tidy_aga_week.first=tidy(mblm(week.first~Year,data=Aga))
#date of last emergence trends
mk_aga_week.last=tidy(mk.test(Aga$week.last))
tidy_aga_week.last=tidy(mblm(week.last~Year,data=Aga))
#duration trends
mk_aga_duration=tidy(mk.test(Aga$duration))
tidy_aga_duration=tidy(mblm(duration~Year,data=Aga))
#multiyear abundance trends
mk_aga_value=tidy(mk.test(Aga$value))
tidy_aga_value=tidy(mblm(value~Year,data=Aga))

##########################################
#same structure for the rest of the taxa


mk_Amp_ct=tidy(mk.test(Amp$ct.week))
tidy_Amp_ct=tidy(mblm(ct.week~Year,data=Amp))


mk_Amp_week.first=tidy(mk.test(Amp$week.first))
tidy_Amp_week.first=tidy(mblm(week.first~Year,data=Amp))


mk_Amp_week.last=tidy(mk.test(Amp$week.last))
tidy_Amp_week.last=tidy(mblm(week.last~Year,data=Amp))


mk_Amp_duration=tidy(mk.test(Amp$duration))
tidy_Amp_duration=tidy(mblm(duration~Year,data=Amp))


mk_Amp_value=tidy(mk.test(Amp$value))
tidy_Amp_value=tidy(mblm(value~Year,data=Amp))

##############################################


mk_Apa_ct=tidy(mk.test(Apa$ct.week))
tidy_Apa_ct=tidy(mblm(ct.week~Year,data=Apa))

mk_Apa_week.first=tidy(mk.test(Apa$week.first))
tidy_Apa_week.first=tidy(mblm(week.first~Year,data=Apa))


mk_Apa_week.last=tidy(mk.test(Apa$week.last))
tidy_Apa_week.last=tidy(mblm(week.last~Year,data=Apa))

mk_Apa_duration=tidy(mk.test(Apa$duration))
tidy_Apa_duration=tidy(mblm(duration~Year,data=Apa))

mk_Apa_value=tidy(mk.test(Apa$value))
tidy_Apa_value=tidy(mblm(value~Year,data=Apa))

###########################################


mk_Bro_ct=tidy(mk.test(Bro$ct.week))
tidy_Bro_ct=tidy(mblm(ct.week~Year,data=Bro))

mk_Bro_week.first=tidy(mk.test(Bro$week.first))
tidy_Bro_week.first=tidy(mblm(week.first~Year,data=Bro))


mk_Bro_week.last=tidy(mk.test(Bro$week.last))
tidy_Bro_week.last=tidy(mblm(week.last~Year,data=Bro))


mk_Bro_duration=tidy(mk.test(Bro$duration))
tidy_Bro_duration=tidy(mblm(duration~Year,data=Bro))

mk_Bro_value=tidy(mk.test(Bro$value))
tidy_Bro_value=tidy(mblm(value~Year,data=Bro))
#mk_Bro_value.boot=mmkh(Bro$value)
#######################################

mk_Bve_ct=tidy(mk.test(Bve$ct.week))
tidy_Bve_ct=tidy(mblm(ct.week~Year,data=Bve))

mk_Bve_week.first=tidy(mk.test(Bve$week.first))
tidy_Bve_week.first=tidy(mblm(week.first~Year,data=Bve))

mk_Bve_week.last=tidy(mk.test(Bve$week.last))
tidy_Bve_week.last=tidy(mblm(week.last~Year,data=Bve))

mk_Bve_duration=tidy(mk.test(Bve$duration))
tidy_Bve_duration=tidy(mblm(duration~Year,data=Bve))


mk_Bve_value=tidy(mk.test(Bve$value))
tidy_Bve_value=tidy(mblm(value~Year,data=Bve))

##########################################


mk_Bra_ct=tidy(mk.test(Bra$ct.week))
tidy_Bra_ct=tidy(mblm(ct.week~Year,data=Bra))

mk_Bra_week.first=tidy(mk.test(Bra$week.first))
tidy_Bra_week.first=tidy(mblm(week.first~Year,data=Bra))

mk_Bra_week.last=tidy(mk.test(Bra$week.last))
tidy_Bra_week.last=tidy(mblm(week.last~Year,data=Bra))

mk_Bra_duration=tidy(mk.test(Bra$duration))
tidy_Bra_duration=tidy(mblm(duration~Year,data=Bra))

mk_Bra_value=tidy(mk.test(Bra$value))
tidy_Bra_value=tidy(mblm(value~Year,data=Bra))
##############################################

mk_Cen_ct=tidy(mk.test(Cen$ct.week))
tidy_Cen_ct=tidy(mblm(ct.week~Year,data=Cen))
#mk_Cen_ct.week.boot=mmkh(Cen$ct.week)

mk_Cen_week.first=tidy(mk.test(Cen$week.first))
tidy_Cen_week.first=tidy(mblm(week.first~Year,data=Cen))
#mk_Cen_week.first_boot=mmkh(Cen$week.first)

mk_Cen_week.last=tidy(mk.test(Cen$week.last))
tidy_Cen_week.last=tidy(mblm(week.last~Year,data=Cen))

mk_Cen_duration=tidy(mk.test(Cen$duration))
tidy_Cen_duration=tidy(mblm(duration~Year,data=Cen))
#mk_Cen_duration_boot=mmkh(Cen$duration)

mk_Cen_value=tidy(mk.test(Cen$value))
tidy_Cen_value=tidy(mblm(value~Year,data=Cen))
#mk_Cen_value_boot=mmkh(Cen$value)
#############################################

mk_Cha_ct=tidy(mk.test(Cha$ct.week))
tidy_Cha_ct=tidy(mblm(ct.week~Year,data=Cha))


mk_Cha_week.first=tidy(mk.test(Cha$week.first))
tidy_Cha_week.first=tidy(mblm(week.first~Year,data=Cha))

mk_Cha_week.last=tidy(mk.test(Cha$week.last))
tidy_Cha_week.last=tidy(mblm(week.last~Year,data=Cha))

mk_Cha_duration=tidy(mk.test(Cha$duration))
tidy_Cha_duration=tidy(mblm(duration~Year,data=Cha))

mk_Cha_value=tidy(mk.test(Cha$value))
tidy_Cha_value=tidy(mblm(value~Year,data=Cha))

##################################################
mk_Dru_ct=tidy(mk.test(Dru$ct.week))
tidy_Dru_ct=tidy(mblm(ct.week~Year,data=Dru))

mk_Dru_week.first=tidy(mk.test(Dru$week.first))
tidy_Dru_week.first=tidy(mblm(week.first~Year,data=Dru))

mk_Dru_week.last=tidy(mk.test(Dru$week.last))
tidy_Dru_week.last=tidy(mblm(week.last~Year,data=Dru))

mk_Dru_duration=tidy(mk.test(Dru$duration))
tidy_Dru_duration=tidy(mblm(duration~Year,data=Dru))


mk_Dru_value=tidy(mk.test(Dru$value))
tidy_Dru_value=tidy(mblm(value~Year,data=Dru))

###############################################



mk_Epi_ct=tidy(mk.test(Epi$ct.week))
tidy_Epi_ct=tidy(mblm(ct.week~Year,data=Epi))

mk_Epi_week.first=tidy(mk.test(Epi$week.first))
tidy_Epi_week.first=tidy(mblm(week.first~Year,data=Epi))

mk_Epi_week.last=tidy(mk.test(Epi$week.last))
tidy_Epi_week.last=tidy(mblm(week.last~Year,data=Epi))


mk_Epi_duration=tidy(mk.test(Epi$duration))
tidy_Epi_duration=tidy(mblm(duration~Year,data=Epi))
#mk_Epi_duration_boot=mmkh(Epi$duration)

mk_Epi_value=tidy(mk.test(Epi$value))
tidy_Epi_value=tidy(mblm(value~Year,data=Epi))

#################################################

mk_Epm_ct=tidy(mk.test(Epm$ct.week))
tidy_Epm_ct=tidy(mblm(ct.week~Year,data=Epm))

mk_Epm_week.first=tidy(mk.test(Epm$week.first))
tidy_Epm_week.first=tidy(mblm(week.first~Year,data=Epm))

mk_Epm_week.last=tidy(mk.test(Epm$week.last))
tidy_Epm_week.last=tidy(mblm(week.last~Year,data=Epm))

mk_Epm_duration=tidy(mk.test(Epm$duration))
tidy_Epm_duration=tidy(mblm(duration~Year,data=Epm))

mk_Epm_value=tidy(mk.test(Epm$value))
tidy_Epm_value=tidy(mblm(value~Year,data=Epm))
##################################################


mk_Iso_ct=tidy(mk.test(Iso$ct.week))
tidy_Iso_ct=tidy(mblm(ct.week~Year,data=Iso))



mk_Iso_week.first=tidy(mk.test(Iso$week.first))
tidy_Iso_week.first=tidy(mblm(week.first~Year,data=Iso))


mk_Iso_week.last=tidy(mk.test(Iso$week.last))
tidy_Iso_week.last=tidy(mblm(week.last~Year,data=Iso))


mk_Iso_duration=tidy(mk.test(Iso$duration))
tidy_Iso_duration=tidy(mblm(duration~Year,data=Iso))


mk_Iso_value=tidy(mk.test(Iso$value))
tidy_Iso_value=tidy(mblm(value~Year,data=Iso))

########################################


mk_Leu_dig_ct=tidy(mk.test(Leu_dig$ct.week))
tidy_Leu_dig_ct=tidy(mblm(ct.week~Year,data=Leu_dig))


mk_Leu_dig_week.first=tidy(mk.test(Leu_dig$week.first))
tidy_Leu_dig_week.first=tidy(mblm(week.first~Year,data=Leu_dig))

mk_Leu_dig_week.last=tidy(mk.test(Leu_dig$week.last))
tidy_Leu_dig_week.last=tidy(mblm(week.last~Year,data=Leu_dig))

mk_Leu_dig_duration=tidy(mk.test(Leu_dig$duration))
tidy_Leu_dig_duration=tidy(mblm(duration~Year,data=Leu_dig))

mk_Leu_dig_value=tidy(mk.test(Leu_dig$value))
tidy_Leu_dig_value=tidy(mblm(value~Year,data=Leu_dig))

###########################################



mk_Leu_nig_ct=tidy(mk.test(Leu_nig$ct.week))
tidy_Leu_nig_ct=tidy(mblm(ct.week~Year,data=Leu_nig))

mk_Leu_nig_week.first=tidy(mk.test(Leu_nig$week.first))
tidy_Leu_nig_week.first=tidy(mblm(week.first~Year,data=Leu_nig))
#mk_Leu_nig_week.first_boot=mmkh(Leu_nig$week.first)

mk_Leu_nig_week.last=tidy(mk.test(Leu_nig$week.last))
tidy_Leu_nig_week.last=tidy(mblm(week.last~Year,data=Leu_nig))

mk_Leu_nig_duration=tidy(mk.test(Leu_nig$duration))
tidy_Leu_nig_duration=tidy(mblm(duration~Year,data=Leu_nig))

mk_Leu_nig_value=tidy(mk.test(Leu_nig$value))
tidy_Leu_nig_value=tidy(mblm(value~Year,data=Leu_nig))
################################################



mk_Leu_pri_ct=tidy(mk.test(Leu_pri$ct.week))
tidy_Leu_pri_ct=tidy(mblm(ct.week~Year,data=Leu_pri))

mk_Leu_pri_week.first=tidy(mk.test(Leu_pri$week.first))
tidy_Leu_pri_week.first=tidy(mblm(week.first~Year,data=Leu_pri))

mk_Leu_pri_week.last=tidy(mk.test(Leu_pri$week.last))
tidy_Leu_pri_week.last=tidy(mblm(week.last~Year,data=Leu_pri))

mk_Leu_pri_duration=tidy(mk.test(Leu_pri$duration))
tidy_Leu_pri_duration=tidy(mblm(duration~Year,data=Leu_pri))

mk_Leu_pri_value=tidy(mk.test(Leu_pri$value))
tidy_Leu_pri_value=tidy(mblm(value~Year,data=Leu_pri))

###############################################



mk_Nem_cam_ct=tidy(mk.test(Nem_cam$ct.week))
tidy_Nem_cam_ct=tidy(mblm(ct.week~Year,data=Nem_cam))


mk_Nem_cam_week.first=tidy(mk.test(Nem_cam$week.first))
tidy_Nem_cam_week.first=tidy(mblm(week.first~Year,data=Nem_cam))


mk_Nem_cam_week.last=tidy(mk.test(Nem_cam$week.last))
tidy_Nem_cam_week.last=tidy(mblm(week.last~Year,data=Nem_cam))


mk_Nem_cam_duration=tidy(mk.test(Nem_cam$duration))
tidy_Nem_cam_duration=tidy(mblm(duration~Year,data=Nem_cam))
#mk_Nem_cam_duration_boot= mmkh(Nem_cam$duration)

mk_Nem_cam_value=tidy(mk.test(Nem_cam$value))
tidy_Nem_cam_value=tidy(mblm(value~Year,data=Nem_cam))

##############################################




mk_Nem_fle_ct=tidy(mk.test(Nem_fle$ct.week))
tidy_Nem_fle_ct=tidy(mblm(ct.week~Year,data=Nem_fle))

mk_Nem_fle_week.first=tidy(mk.test(Nem_fle$week.first))
tidy_Nem_fle_week.first=tidy(mblm(week.first~Year,data=Nem_fle))

mk_Nem_fle_week.last=tidy(mk.test(Nem_fle$week.last))
tidy_Nem_fle_week.last=tidy(mblm(week.last~Year,data=Nem_fle))

mk_Nem_fle_duration=tidy(mk.test(Nem_fle$duration))
tidy_Nem_fle_duration=tidy(mblm(duration~Year,data=Nem_fle))

mk_Nem_fle_value=tidy(mk.test(Nem_fle$value))
tidy_Nem_fle_value=tidy(mblm(value~Year,data=Nem_fle))
###############################################


mk_Nem_mar_ct=tidy(mk.test(Nem_mar$ct.week))
tidy_Nem_mar_ct=tidy(mblm(ct.week~Year,data=Nem_mar))

mk_Nem_mar_week.first=tidy(mk.test(Nem_mar$week.first))
tidy_Nem_mar_week.first=tidy(mblm(week.first~Year,data=Nem_mar))

mk_Nem_mar_week.last=tidy(mk.test(Nem_mar$week.last))
tidy_Nem_mar_week.last=tidy(mblm(week.last~Year,data=Nem_mar))

mk_Nem_mar_duration=tidy(mk.test(Nem_mar$duration))
tidy_Nem_mar_duration=tidy(mblm(duration~Year,data=Nem_mar))


mk_Nem_mar_value=tidy(mk.test(Nem_mar$value))
tidy_Nem_mar_value=tidy(mblm(value~Year,data=Nem_mar))

################################################

mk_Nrl_pic_ct=tidy(mk.test(Nrl_pic$ct.week))
tidy_Nrl_pic_ct=tidy(mblm(ct.week~Year,data=Nrl_pic))


mk_Nrl_pic_week.first=tidy(mk.test(Nrl_pic$week.first))
tidy_Nrl_pic_week.first=tidy(mblm(week.first~Year,data=Nrl_pic))

mk_Nrl_pic_week.last=tidy(mk.test(Nrl_pic$week.last))
tidy_Nrl_pic_week.last=tidy(mblm(week.last~Year,data=Nrl_pic))


mk_Nrl_pic_duration=tidy(mk.test(Nrl_pic$duration))
tidy_Nrl_pic_duration=tidy(mblm(duration~Year,data=Nrl_pic))


mk_Nrl_pic_value=tidy(mk.test(Nrl_pic$value))
tidy_Nrl_pic_value=tidy(mblm(value~Year,data=Nrl_pic))

################################################
mk_Par_sub_ct=tidy(mk.test(Par_sub$ct.week))
tidy_Par_sub_ct=tidy(mblm(ct.week~Year,data=Par_sub))

mk_Par_sub_week.first=tidy(mk.test(Par_sub$week.first))
tidy_Par_sub_week.first=tidy(mblm(week.first~Year,data=Par_sub))


mk_Par_sub_week.last=tidy(mk.test(Par_sub$week.last))
tidy_Par_sub_week.last=tidy(mblm(week.last~Year,data=Par_sub))

mk_Par_sub_duration=tidy(mk.test(Par_sub$duration))
tidy_Par_sub_duration=tidy(mblm(duration~Year,data=Par_sub))

mk_Par_sub_value=tidy(mk.test(Par_sub$value))
tidy_Par_sub_value=tidy(mblm(value~Year,data=Par_sub))

#########################################

mk_Ple_con_ct=tidy(mk.test(Ple_con$ct.week))
tidy_Ple_con_ct=tidy(mblm(ct.week~Year,data=Ple_con))

mk_Ple_con_week.first=tidy(mk.test(Ple_con$week.first))
tidy_Ple_con_week.first=tidy(mblm(week.first~Year,data=Ple_con))

mk_Ple_con_week.last=tidy(mk.test(Ple_con$week.last))
tidy_Ple_con_week.last=tidy(mblm(week.last~Year,data=Ple_con))

mk_Ple_con_duration=tidy(mk.test(Ple_con$duration))
tidy_Ple_con_duration=tidy(mblm(duration~Year,data=Ple_con))

mk_Ple_con_value=tidy(mk.test(Ple_con$value))
tidy_Ple_con_value=tidy(mblm(value~Year,data=Ple_con))

#########################################################

mk_Pot_luc_ct=tidy(mk.test(Pot_luc$ct.week))
tidy_Pot_luc_ct=tidy(mblm(ct.week~Year,data=Pot_luc))

mk_Pot_luc_week.first=tidy(mk.test(Pot_luc$week.first))
tidy_Pot_luc_week.first=tidy(mblm(week.first~Year,data=Pot_luc))
#mk_Pot_luc_week.first.boot=mmkh(Pot_luc$week.first)

mk_Pot_luc_week.last=tidy(mk.test(Pot_luc$week.last))
tidy_Pot_luc_week.last=tidy(mblm(week.last~Year,data=Pot_luc))

mk_Pot_luc_duration=tidy(mk.test(Pot_luc$duration))
tidy_Pot_luc_duration=tidy(mblm(duration~Year,data=Pot_luc))

mk_Pot_luc_value=tidy(mk.test(Pot_luc$value))
tidy_Pot_luc_value=tidy(mblm(value~Year,data=Pot_luc))

####################################################


mk_Pro_aub_ct=tidy(mk.test(Pro_aub$ct.week))
tidy_Pro_aub_ct=tidy(mblm(ct.week~Year,data=Pro_aub))


mk_Pro_aub_week.first=tidy(mk.test(Pro_aub$week.first))
tidy_Pro_aub_week.first=tidy(mblm(week.first~Year,data=Pro_aub))


mk_Pro_aub_week.last=tidy(mk.test(Pro_aub$week.last))
tidy_Pro_aub_week.last=tidy(mblm(week.last~Year,data=Pro_aub))
#mk_Pro_aub_week.last_boot=mmkh(Pro_aub$week.last)


mk_Pro_aub_duration=tidy(mk.test(Pro_aub$duration))
tidy_Pro_aub_duration=tidy(mblm(duration~Year,data=Pro_aub))

mk_Pro_aub_value=tidy(mk.test(Pro_aub$value))
tidy_Pro_aub_value=tidy(mblm(value~Year,data=Pro_aub))

#############################################


mk_Pro_int_ct=tidy(mk.test(Pro_int$ct.week))
tidy_Pro_int_ct=tidy(mblm(ct.week~Year,data=Pro_int))

mk_Pro_int_week.first=tidy(mk.test(Pro_int$week.first))
tidy_Pro_int_week.first=tidy(mblm(week.first~Year,data=Pro_int))

mk_Pro_int_week.last=tidy(mk.test(Pro_int$week.last))
tidy_Pro_int_week.last=tidy(mblm(week.last~Year,data=Pro_int))

mk_Pro_int_duration=tidy(mk.test(Pro_int$duration))
tidy_Pro_int_duration=tidy(mblm(duration~Year,data=Pro_int))

mk_Pro_int_value=tidy(mk.test(Pro_int$value))
tidy_Pro_int_value=tidy(mblm(value~Year,data=Pro_int))
#mk_Pro_int_value.boot=mmkh(Pro_int$value)
###################################################


mk_Pro_mey_ct=tidy(mk.test(Pro_mey$ct.week))
tidy_Pro_mey_ct=tidy(mblm(ct.week~Year,data=Pro_mey))

mk_Pro_mey_week.first=tidy(mk.test(Pro_mey$week.first))
tidy_Pro_mey_week.first=tidy(mblm(week.first~Year,data=Pro_mey))

mk_Pro_mey_week.last=tidy(mk.test(Pro_mey$week.last))
tidy_Pro_mey_week.last=tidy(mblm(week.last~Year,data=Pro_mey))

mk_Pro_mey_duration=tidy(mk.test(Pro_mey$duration))
tidy_Pro_mey_duration=tidy(mblm(duration~Year,data=Pro_mey))

mk_Pro_mey_value=tidy(mk.test(Pro_mey$value))
tidy_Pro_mey_value=tidy(mblm(value~Year,data=Pro_mey))
###########################################################
mk_Rhy_fas_ct=tidy(mk.test(Rhy_fas$ct.week))
tidy_Rhy_fas_ct=tidy(mblm(ct.week~Year,data=Rhy_fas))

mk_Rhy_fas_week.first=tidy(mk.test(Rhy_fas$week.first))
tidy_Rhy_fas_week.first=tidy(mblm(week.first~Year,data=Rhy_fas))

mk_Rhy_fas_week.last=tidy(mk.test(Rhy_fas$week.last))
tidy_Rhy_fas_week.last=tidy(mblm(week.last~Year,data=Rhy_fas))

mk_Rhy_fas_duration=tidy(mk.test(Rhy_fas$duration))
tidy_Rhy_fas_duration=tidy(mblm(duration~Year,data=Rhy_fas))

mk_Rhy_fas_value=tidy(mk.test(Rhy_fas$value))
tidy_Rhy_fas_value=tidy(mblm(value~Year,data=Rhy_fas))

################################################

mk_Ser_per_ct=tidy(mk.test(Ser_per$ct.week))
tidy_Ser_per_ct=tidy(mblm(ct.week~Year,data=Ser_per))

mk_Ser_per_week.first=tidy(mk.test(Ser_per$week.first))
tidy_Ser_per_week.first=tidy(mblm(week.first~Year,data=Ser_per))

mk_Ser_per_week.last=tidy(mk.test(Ser_per$week.last))
tidy_Ser_per_week.last=tidy(mblm(week.last~Year,data=Ser_per))

mk_Ser_per_duration=tidy(mk.test(Ser_per$duration))
tidy_Ser_per_duration=tidy(mblm(duration~Year,data=Ser_per))

mk_Ser_per_value=tidy(mk.test(Ser_per$value))
tidy_Ser_per_value=tidy(mblm(value~Year,data=Ser_per))
#mk_Ser_per_value.boot=mmkh(Ser_per$value)
################################################

mk_Sil_pal_ct=tidy(mk.test(Sil_pal$ct.week))
tidy_Sil_pal_ct=tidy(mblm(ct.week~Year,data=Sil_pal))

mk_Sil_pal_week.first=tidy(mk.test(Sil_pal$week.first))
tidy_Sil_pal_week.first=tidy(mblm(week.first~Year,data=Sil_pal))

mk_Sil_pal_week.last=tidy(mk.test(Sil_pal$week.last))
tidy_Sil_pal_week.last=tidy(mblm(week.last~Year,data=Sil_pal))

mk_Sil_pal_duration=tidy(mk.test(Sil_pal$duration))
tidy_Sil_pal_duration=tidy(mblm(duration~Year,data=Sil_pal))


mk_Sil_pal_value=tidy(mk.test(Sil_pal$value))
tidy_Sil_pal_value=tidy(mblm(value~Year,data=Sil_pal))

########################################################


mk_Sip_tor_ct=tidy(mk.test(Sip_tor$ct.week))
tidy_Sip_tor_ct=tidy(mblm(ct.week~Year,data=Sip_tor))

mk_Sip_tor_week.first=tidy(mk.test(Sip_tor$week.first))
tidy_Sip_tor_week.first=tidy(mblm(week.first~Year,data=Sip_tor))

mk_Sip_tor_week.last=tidy(mk.test(Sip_tor$week.last))
tidy_Sip_tor_week.last=tidy(mblm(week.last~Year,data=Sip_tor))
mk_Sip_tor_week.last.boot=tsboot(Sip_tor$week.last, MKtau, R=500, l=5, sim="fixed")

mk_Sip_tor_duration=tidy(mk.test(Sip_tor$duration))
tidy_Sip_tor_duration=tidy(mblm(duration~Year,data=Sip_tor))

mk_Sip_tor_value=tidy(mk.test(Sip_tor$value))
tidy_Sip_tor_value=tidy(mblm(value~Year,data=Sip_tor))
#mk_Sip_tor_value.boot=mmkh(Sip_tor$value)
##################################################


mk_Tin_ros_ct=tidy(mk.test(Tin_ros$ct.week))
tidy_Tin_ros_ct=tidy(mblm(ct.week~Year,data=Tin_ros))

mk_Tin_ros_week.first=tidy(mk.test(Tin_ros$week.first))
tidy_Tin_ros_week.first=tidy(mblm(week.first~Year,data=Tin_ros))

mk_Tin_ros_week.last=tidy(mk.test(Tin_ros$week.last))
tidy_Tin_ros_week.last=tidy(mblm(week.last~Year,data=Tin_ros))

mk_Tin_ros_duration=tidy(mk.test(Tin_ros$duration))
tidy_Tin_ros_duration=tidy(mblm(duration~Year,data=Tin_ros))

mk_Tin_ros_value=tidy(mk.test(Tin_ros$value))
tidy_Tin_ros_value=tidy(mblm(value~Year,data=Tin_ros))

##################################################

mk_Wor_occ_ct=tidy(mk.test(Wor_occ$ct.week))
tidy_Wor_occ_ct=tidy(mblm(ct.week~Year,data=Wor_occ))

mk_Wor_occ_week.first=tidy(mk.test(Wor_occ$week.first))
tidy_Wor_occ_week.first=tidy(mblm(week.first~Year,data=Wor_occ))


mk_Wor_occ_week.last=tidy(mk.test(Wor_occ$week.last))
tidy_Wor_occ_week.last=tidy(mblm(week.last~Year,data=Wor_occ))

mk_Wor_occ_duration=tidy(mk.test(Wor_occ$duration))
tidy_Wor_occ_duration=tidy(mblm(duration~Year,data=Wor_occ))
#mk_Wor_occ_duration.boot=mmkh(Wor_occ$duration)

mk_Wor_occ_value=tidy(mk.test(Wor_occ$value))
tidy_Wor_occ_value=tidy(mblm(value~Year,data=Wor_occ))



########################################################

#extract Sen's slope and Mann Kendall-test coefficients as dataframes

#sens slopes combined
tidyy=bind_rows(tidy_aga_ct	,
                tidy_aga_duration	,
                tidy_aga_week.last	,
                tidy_aga_week.first	,
                tidy_aga_value,
                tidy_Amp_ct	,
                tidy_Amp_duration	,
                tidy_Amp_week.first	,
                tidy_Amp_week.last	,
                tidy_Amp_value,
                tidy_Apa_ct	,
                tidy_Apa_duration	,
                tidy_Apa_week.first	,
                tidy_Apa_week.last	,
                tidy_Apa_value,
                tidy_Bra_ct	,
                tidy_Bra_duration	,
                tidy_Bra_week.first	,
                tidy_Bra_week.last	,
                tidy_Bra_value,
                tidy_Bro_ct	,
                tidy_Bro_duration	,
                tidy_Bro_week.first	,
                tidy_Bro_week.last	,
                tidy_Bro_value,
                tidy_Bve_ct	,
                tidy_Bve_duration	,
                tidy_Bve_week.first	,
                tidy_Bve_week.last	,
                tidy_Bve_value,
                tidy_Cen_ct	,
                tidy_Cen_duration	,
                tidy_Cen_week.first	,
                tidy_Cen_week.last	,
                tidy_Cen_value	,
                tidy_Cha_ct	,
                tidy_Cha_duration	,
                tidy_Cha_week.first	,
                tidy_Cha_week.last	,
                tidy_Cha_value,
                tidy_Dru_ct	,
                tidy_Dru_duration	,
                tidy_Dru_week.first	,
                tidy_Dru_week.last	,
                tidy_Dru_value,
                tidy_Epi_ct	,
                tidy_Epi_duration	,
                tidy_Epi_week.first	,
                tidy_Epi_week.last	,
                tidy_Epi_value,
                tidy_Epm_ct	,
                tidy_Epm_duration	,
                tidy_Epm_week.first	,
                tidy_Epm_week.last	,
                tidy_Epm_value,
                tidy_Iso_ct	,
                tidy_Iso_duration	,
                tidy_Iso_week.first	,
                tidy_Iso_week.last	,
                tidy_Iso_value,
                tidy_Leu_dig_ct	,
                tidy_Leu_dig_duration	,
                tidy_Leu_dig_week.first	,
                tidy_Leu_dig_week.last	,
                tidy_Leu_dig_value,
                tidy_Leu_nig_ct	,
                tidy_Leu_nig_duration	,
                tidy_Leu_nig_week.first	,
                tidy_Leu_nig_week.last	,
                tidy_Leu_nig_value,
                tidy_Leu_pri_ct	,
                tidy_Leu_pri_duration	,
                tidy_Leu_pri_week.first	,
                tidy_Leu_pri_week.last	,
                tidy_Leu_pri_value,
                tidy_Nem_cam_ct	,
                tidy_Nem_cam_duration	,
                tidy_Nem_cam_week.first	,
                tidy_Nem_cam_week.last	,
                tidy_Nem_cam_value,
                tidy_Nem_fle_ct	,
                tidy_Nem_fle_duration	,
                tidy_Nem_fle_week.first	,
                tidy_Nem_fle_week.last	,
                tidy_Nem_fle_value,
                tidy_Nem_mar_ct	,
                tidy_Nem_mar_duration	,
                tidy_Nem_mar_week.first	,
                tidy_Nem_mar_week.last	,
                tidy_Nem_mar_value,
                tidy_Nrl_pic_ct	,
                tidy_Nrl_pic_duration	,
                tidy_Nrl_pic_week.first	,
                tidy_Nrl_pic_week.last	,
                tidy_Nrl_pic_value,
                tidy_Par_sub_ct	,
                tidy_Par_sub_duration	,
                tidy_Par_sub_week.first	,
                tidy_Par_sub_week.last	,
                tidy_Par_sub_value,
                tidy_Ple_con_ct	,
                tidy_Ple_con_duration	,
                tidy_Ple_con_week.first	,
                tidy_Ple_con_week.last	,
                tidy_Ple_con_value,
                tidy_Pot_luc_ct	,
                tidy_Pot_luc_duration	,
                tidy_Pot_luc_week.first	,
                tidy_Pot_luc_week.last	,
                tidy_Pot_luc_value,
                tidy_Pro_aub_ct	,
                tidy_Pro_aub_duration	,
                tidy_Pro_aub_week.first	,
                tidy_Pro_aub_week.last	,
                tidy_Pro_aub_value,
                tidy_Pro_int_ct	,
                tidy_Pro_int_duration	,
                tidy_Pro_int_week.first	,
                tidy_Pro_int_week.last	,
                tidy_Pro_int_value,
                tidy_Pro_mey_ct	,
                tidy_Pro_mey_duration	,
                tidy_Pro_mey_week.first	,
                tidy_Pro_mey_week.last	,
                tidy_Pro_mey_value,
                tidy_Rhy_fas_ct	,
                tidy_Rhy_fas_duration	,
                tidy_Rhy_fas_week.first	,
                tidy_Rhy_fas_week.last	,
                tidy_Rhy_fas_value,
                tidy_Ser_per_ct	,
                tidy_Ser_per_duration	,
                tidy_Ser_per_week.first	,
                tidy_Ser_per_week.last	,
                tidy_Ser_per_value,
                tidy_Sil_pal_ct	,
                tidy_Sil_pal_duration	,
                tidy_Sil_pal_week.first	,
                tidy_Sil_pal_week.last	,
                tidy_Sil_pal_value,
                tidy_Sip_tor_ct	,
                tidy_Sip_tor_duration	,
                tidy_Sip_tor_week.first	,
                tidy_Sip_tor_week.last	,
                tidy_Sip_tor_value,
                tidy_Tin_ros_ct	,
                tidy_Tin_ros_duration	,
                tidy_Tin_ros_week.first	,
                tidy_Tin_ros_week.last,
                tidy_Tin_ros_value,
                tidy_Wor_occ_ct	,
                tidy_Wor_occ_duration	,
                tidy_Wor_occ_week.first	,
                tidy_Wor_occ_week.last,
                tidy_Wor_occ_value)

#add names of the taxa and variable to the dataset

name=c("aga.ct",
       "aga.duration",
       "aga.week.first",
       "aga.week.last",
       "aga.value",
       "Amp.ct",
       "Amp.duration",
       "Amp.week.first",
       "Amp.week.last",
       "Amp.value",
       "Apa.ct",
       "Apa.duration",
       "Apa.week.first",
       "Apa.week.last",
       "Apa.value",
       "Bra.ct",
       "Bra.duration",
       "Bra.week.first",
       "Bra.week.last",
       "Bra.value",
       "Bro.ct",
       "Bro.duration",
       "Bro.week.first",
       "Bro.week.last",
       "Bro.value",
       "Bve.ct",
       "Bve.duration",
       "Bve.week.first",
       "Bve.week.last",
       "Bve.value",
       "Cen.ct",
       "Cen.duration",
       "Cen.week.first",
       "Cen.week.last",
       "Cen.value",
       "Cha.ct",
       "Cha.duration",
       "Cha.week.first",
       "Cha.week.last",
       "Cha.value",
       "Dru.ct",
       "Dru.duration",
       "Dru.week.first",
       "Dru.week.last",
       "Dru.value",
       "Epi.ct",
       "Epi.duration",
       "Epi.week.first",
       "Epi.week.last",
       "Epi.value",
       "Epm.ct",
       "Epm.duration",
       "Epm.week.first",
       "Epm.week.last",
       "Epm.value",
       "Iso.ct",
       "Iso.duration",
       "Iso.week.first",
       "Iso.week.last",
       "Iso.value",
       "Leu.dig.ct",
       "Leu.dig.duration",
       "Leu.dig.first",
       "Leu.dig.last",
       "Leu.dig.value",
       "Leu.nig.ct",
       "Leu.nig.duration",
       "Leu.nig.first",
       "Leu.nig.last",
       "Leu.nig.value",
       "Leu.pri.ct",
       "Leu.pri.duration",
       "Leu.pri.first",
       "Leu.pri.last",
       "Leu.pri.value",
       "Nem.cam.ct",
       "Nem.cam.duration",
       "Nem.cam.first",
       "Nem.cam.last",
       "Nem.cam.value",
       "Nem.fle.ct",
       "Nem.fle.duration",
       "Nem.fle.first",
       "Nem.fle.last",
       "Nem.fle.value",
       "Nem.mar.ct",
       "Nem.mar.duration",
       "Nem.mar.first",
       "Nem.mar.last",
       "Nem.mar.value",
       "Nrl.pic.ct",
       "Nrl.pic.duration",
       "Nrl.pic.first",
       "Nrl.pic.last",
       "Nrl.pic.value",
       "Par.sub.ct",
       "Par.sub.duration",
       "Par.sub.first",
       "Par.sub.last",
       "Par.sub.value",
       "Ple.con.ct",
       "Ple.con.duration",
       "Ple.con.first",
       "Ple.con.last",
       "Ple.con.value",
       "Pot.luc.ct",
       "Pot.luc.duration",
       "Pot.luc.first",
       "Pot.luc.last",
       "Pot.luc.value",
       "Pro.aub.ct",
       "Pro.aub.duration",
       "Pro.aub.first",
       "Pro.aub.last",
       "Pro.aub.value",
       "Pro.int.ct",
       "Pro.int.duration",
       "Pro.int.first",
       "Pro.int.last",
       "Pro.int.value",
       "Pro.mey.ct",
       "Pro.mey.duration",
       "Pro.mey.first",
       "Pro.mey.last",
       "Pro.mey.value",
       "Rhy.fas.ct",
       "Rhy.fas.duration",
       "Rhy.fas.first",
       "Rhy.fas.last",
       "Rhy.fas.value",
       "Ser.per.ct",
       "Ser.per.duration",
       "Ser.per.first",
       "Ser.per.last",
       "Ser.per.value",
       "Sil.pal.ct",
       "Sil.pal.duration",
       "Sil.pal.first",
       "Sil.pal.last",
       "Sil.pal.value",
       "Sip.tor.ct",
       "Sip.tor.duration",
       "Sip.tor.first",
       "Sip.tor.last",
       "Sip.tor.value",
       "Tin.ros.ct",
       "Tin.ros.duration",
       "Tin.ros.first",
       "Tin.ros.last",
       "Tin.ros.value",
       "Wor.occ.ct",
       "Wor.occ.duration",
       "Wor.occ.first",
       "Wor.occ.last",
       "Wor.occ.value")


names=rep(name,each=2) #creating vector of names for the coefs dataframe

df1=cbind(names,tidyy)

#same for the MK-test statistics

tidmk=bind_rows(mk_aga_ct	,
                mk_aga_duration	,
                mk_aga_week.first	,
                mk_aga_week.last	,
                mk_aga_value,
                mk_Amp_ct	,
                mk_Amp_duration	,
                mk_Amp_week.first	,
                mk_Amp_week.last	,
                mk_Amp_value,
                mk_Apa_ct	,
                mk_Apa_duration	,
                mk_Apa_week.first	,
                mk_Apa_week.last	,
                mk_Apa_value,
                mk_Bra_ct	,
                mk_Bra_duration	,
                mk_Bra_week.first	,
                mk_Bra_week.last	,
                mk_Bra_value,
                mk_Bro_ct	,
                mk_Bro_duration	,
                mk_Bro_week.first	,
                mk_Bro_week.last	,
                mk_Bro_value,
                mk_Bve_ct	,
                mk_Bve_duration	,
                mk_Bve_week.first	,
                mk_Bve_week.last	,
                mk_Bve_value,
                mk_Cen_ct	,
                mk_Cen_duration	,
                mk_Cen_week.first	,
                mk_Cen_week.last	,
                mk_Cen_value	,
                mk_Cha_ct	,
                mk_Cha_duration	,
                mk_Cha_week.first	,
                mk_Cha_week.last	,
                mk_Cha_value,
                mk_Dru_ct	,
                mk_Dru_duration	,
                mk_Dru_week.first	,
                mk_Dru_week.last	,
                mk_Dru_value,
                mk_Epi_ct	,
                mk_Epi_duration	,
                mk_Epi_week.first	,
                mk_Epi_week.last	,
                mk_Epi_value,
                mk_Epm_ct	,
                mk_Epm_duration	,
                mk_Epm_week.first	,
                mk_Epm_week.last	,
                mk_Epm_value,
                mk_Iso_ct	,
                mk_Iso_duration	,
                mk_Iso_week.first	,
                mk_Iso_week.last	,
                mk_Iso_value,
                mk_Leu_dig_ct	,
                mk_Leu_dig_duration	,
                mk_Leu_dig_week.first	,
                mk_Leu_dig_week.last	,
                mk_Leu_dig_value,
                mk_Leu_nig_ct	,
                mk_Leu_nig_duration	,
                mk_Leu_nig_week.first	,
                mk_Leu_nig_week.last	,
                mk_Leu_nig_value,
                mk_Leu_pri_ct	,
                mk_Leu_pri_duration	,
                mk_Leu_pri_week.first	,
                mk_Leu_pri_week.last	,
                mk_Leu_pri_value,
                mk_Nem_cam_ct	,
                mk_Nem_cam_duration	,
                mk_Nem_cam_week.first	,
                mk_Nem_cam_week.last	,
                mk_Nem_cam_value,
                mk_Nem_fle_ct	,
                mk_Nem_fle_duration	,
                mk_Nem_fle_week.first	,
                mk_Nem_fle_week.last	,
                mk_Nem_fle_value,
                mk_Nem_mar_ct	,
                mk_Nem_mar_duration	,
                mk_Nem_mar_week.first	,
                mk_Nem_mar_week.last	,
                mk_Nem_mar_value,
                mk_Nrl_pic_ct	,
                mk_Nrl_pic_duration	,
                mk_Nrl_pic_week.first	,
                mk_Nrl_pic_week.last	,
                mk_Nrl_pic_value,
                mk_Par_sub_ct	,
                mk_Par_sub_duration	,
                mk_Par_sub_week.first	,
                mk_Par_sub_week.last	,
                mk_Par_sub_value,
                mk_Ple_con_ct	,
                mk_Ple_con_duration	,
                mk_Ple_con_week.first	,
                mk_Ple_con_week.last	,
                mk_Ple_con_value,
                mk_Pot_luc_ct	,
                mk_Pot_luc_duration	,
                mk_Pot_luc_week.first	,
                mk_Pot_luc_week.last	,
                mk_Pot_luc_value,
                mk_Pro_aub_ct	,
                mk_Pro_aub_duration	,
                mk_Pro_aub_week.first	,
                mk_Pro_aub_week.last	,
                mk_Pro_aub_value,
                mk_Pro_int_ct	,
                mk_Pro_int_duration	,
                mk_Pro_int_week.first	,
                mk_Pro_int_week.last	,
                mk_Pro_int_value,
                mk_Pro_mey_ct	,
                mk_Pro_mey_duration	,
                mk_Pro_mey_week.first	,
                mk_Pro_mey_week.last	,
                mk_Pro_mey_value,
                mk_Rhy_fas_ct	,
                mk_Rhy_fas_duration	,
                mk_Rhy_fas_week.first	,
                mk_Rhy_fas_week.last	,
                mk_Rhy_fas_value,
                mk_Ser_per_ct	,
                mk_Ser_per_duration	,
                mk_Ser_per_week.first	,
                mk_Ser_per_week.last	,
                mk_Ser_per_value,
                mk_Sil_pal_ct	,
                mk_Sil_pal_duration	,
                mk_Sil_pal_week.first	,
                mk_Sil_pal_week.last	,
                mk_Sil_pal_value,
                mk_Sip_tor_ct	,
                mk_Sip_tor_duration	,
                mk_Sip_tor_week.first	,
                mk_Sip_tor_week.last	,
                mk_Sip_tor_value,
                mk_Tin_ros_ct	,
                mk_Tin_ros_duration	,
                mk_Tin_ros_week.first	,
                mk_Tin_ros_week.last,
                mk_Tin_ros_value,
                mk_Wor_occ_ct	,
                mk_Wor_occ_duration	,
                mk_Wor_occ_week.first	,
                mk_Wor_occ_week.last,
                mk_Wor_occ_value)

taxaEPT=c("Aga_fus"	,
          "Amp_sta"	,
          "Apa_fim"	,
          "Bra_ris",
          "Bae_roh"	,
          "Bae_ver",
          "Cen_lut"	,
          "Cha_vil"	,
          "Dru_ann"	,
          "Eph_ign"	,
          "Eph_muc"	,
          "Iso_goe"	,
          "Leu_dig"	,
          "Leu_nig"	,
          "Leu_pri"	,
          "Nem_cam"	,
          "Nem_fle"	,
          "Nem_mar"	,
          "Nrl_pic"	,
          "Par_sub"	,
          "Ple_con"	,
          "Pot_luc"	,
          "Pro_aub"	,
          "Pro_int"	,
          "Pro_mey"	,
          "Rhy_fas"	,
          "Ser_per"	,
          "Sil_pal"	,
          "Sip_tor"	,
          "Tin_ros"	,
          "Wor_occ")

taxa1EPT=c("Agapetus fuscipies"	,
           "Amphinemura standfussi"	,
           "Apatania fimbriata"	,
           "Brachyptera risi",
           "Baetis rohdani"	,
           "Baetis vernus",
           "Centroptilum luteolum"	,
           "Chaetopteryx villosa"	,
           "Drusus annulatus"	,
           "Ephemerella ignita"	,
           "Ephemerella mucronata"	,
           "Isoperla goertzi"	,
           "Leuctra digitata"	,
           "Leuctra nigra"	,
           "Leuctra prima"	,
           "Nemoura cambrica"	,
           "Nemoura flexuosa"	,
           "Nemoura marginata"	,
           "Nemurella  pictetii"	,
           "Paraleptophlebia submarginata"	,
           "Plectrocnemia conspersa"	,
           "Potamophylax luctuosus"	,
           "Protonemura auberti"	,
           "Protonemura intricata"	,
           "Protonemura meyeri"	,
           "Rhyacophila fasciata"	,
           "Sericostoma personatum"	,
           "Silo pallipes"	,
           "Siphonoperla torrentium"	,
           "Tinodes rostocki"	,
           "Wormaldia occipitalis")



#names of taxa for the Sen's slope aggregated
taxas1=rep(taxaEPT,each=5)
#names of taxa for the Mann Kendall aggregated
taxas2=rep(taxaEPT,each=10)

vars=c("ct","duration","first","last","value")#label phenophases
var=rep(vars,31)
#attach taxa names to the Mann Kendall aggregated
df2=cbind(name,tidmk)
tidmk12=tidmk #rename for convinience
#values for autocorrelated time series
#attach the taxa name
df2$spp=taxas1
name=data.frame(name)
#attach taxa names to the sen's slopes dataset
df1=cbind(tidyy,taxas2)
vars1=c("ct","ct","duration","duration","first","first","last","last","value","value")
#attach names of the phenophase
var2=rep(vars1,31)
df1$phenophase=vars1
name$phenophase=var
tidmk12$phenophase=var
#attach names of the phenophase
tidmk12$name=name$name
#attach taxa name

df2<-merge(tidmk12,name, by="name", all.x=TRUE, all.y=TRUE)

#sort significan and insignificant trends

a <- df2 %>%
  filter( p.value>0.05) %>%
  mutate(final = ifelse(p.value >0.05, "insignificant"))
days=rep(0,68)


b <- df2 %>%
  filter( p.value<=0.05) %>%
  mutate(final = ifelse(p.value<=0.05, "significant"))
#re-aggregate

df2=rbind(a,b) # all coefficients

#attach taxa, subset only slopes from Sen's slope
df1$spp=df1$taxas
df1=subset(df1,df1$term=="Year")
df1$name=name$name

#merge Sen's slope and Mann Kendall data
dfx<-merge(df2,df1,by = c("name"), all.x = TRUE, all.y = TRUE)
#sort significan and insignificant trends again

a=subset(dfx,dfx$final=="insignificant")
days=rep(0,87)
a$days=days

b=subset(dfx,dfx$final=="significant")
b$days=b$estimate*7
#subset slopes only

dft=rbind(a,b) # all coefficients
dft$day1=dft$estimate*7 # recalculate number of the day of the 
#year from the number of the week
#remove phenophases not to be displayed
dft1=dft[! dft$phenophase %in% c("last", "first","value"), ]
dft2=dft[! dft$phenophase %in% c("last", "first","ct","duration"), ]
#reordering columns to achive taxa being sorted by the 
#insects orders [Ephemeroptera, trichoptera & Plecoptera]
taxa_com=cbind(taxaEPT,taxa1EPT)
taxa_com=data.frame(taxa_com)
taxa_com$spp=taxa_com$taxaEPT
dft1<-merge(taxa_com,dft1,by = c("spp"), all.x = TRUE, all.y = TRUE)
orders=c("T","T","P","P","T","T","E","E","E","E","P","P","E","E","T","T","T","T","E","E","E","E","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","E","E","T","T","T","T","P","P","P","P","P","P","T","T","T","T","T","T","P","P","T","T","T","T")
dft1$orders=orders
dft1<-with(dft1, dft1[order(orders,taxa1EPT),])

dft2=dft1[,c(3,12,19,21,22)]
dft2$day1<- reorder(dft2$orders,dft2$taxa1EPT)
dft2=dft2[,c(1:3,5)]
dft2$day1=dft1$day1

dft2$taxa1EPT <- factor(dft2$taxa1EPT, levels = c("Baetis rohdani","Baetis vernus","Centroptilum luteolum","Ephemerella ignita","Ephemerella mucronata","Paraleptophlebia submarginata","Amphinemura standfussi","Brachyptera risi","Isoperla goertzi","Leuctra digitata",
                                                  "Leuctra nigra",
                                                  "Leuctra prima",
                                                  "Nemoura cambrica",
                                                  "Nemoura flexuosa",
                                                  "Nemoura marginata",
                                                  "Nemurella  pictetii",           
                                                  "Protonemura auberti",
                                                  "Protonemura intricata",
                                                  "Protonemura meyeri",
                                                  "Siphonoperla torrentium",
                                                  "Agapetus fuscipies",
                                                  "Apatania fimbriata",
                                                  "Chaetopteryx villosa",
                                                  "Drusus annulatus",
                                                  "Plectrocnemia conspersa",
                                                  "Potamophylax luctuosus",
                                                  "Rhyacophila fasciata",
                                                  "Sericostoma personatum",
                                                  "Silo pallipes"	,
                                                  "Tinodes rostocki",           
                                                  "Wormaldia occipitalis"))








g1=ggplot(dft2, aes(as.factor(taxa1EPT), day1))+facet_wrap(~phenophase,scales = "free_x", ncol = 3, switch = "y") + geom_point(size=10,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))+scale_size_manual(values =c("significant" = 6,"insignificant" = 4))+coord_flip()+ geom_hline(yintercept=0,lty=2)+ ylab("[shift, days per year]") +theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))            

#g1=ggplot(dft1, aes(taxa1EPT, day1)) + geom_point(size=10,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))+scale_size_manual(values =c("significant" = 6,"insignificant" = 4))+coord_flip()+ geom_hline(yintercept=0,lty=2)+ ylab("[shift, days per year]") +theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))            
#g1t=g1+facet_wrap(~phenophase,scales = "free_x", ncol = 3, switch = "y") +ylab(NULL) +theme(strip.background = element_blank())+ geom_point(size=5,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))

g1t=g1+theme_bw()
g1t=g1t+theme(legend.position="top")
g1t=g1t+theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))+ theme(text = element_text(size = 40))

italic.10.text <- element_text(face = "italic", size = 30)
#figure5
g1t=g1t+ theme(axis.text.y = italic.10.text)+coord_flip()+theme(strip.background = element_blank(), strip.text = element_blank())



#######################
taxa1EPT=c("Agapetus fuscipies"	,
           "Amphinemura standfussi"	,
           "Apatania fimbriata"	,
           "Brachyptera risi",
           "Baetis rohdani"	,
           "Baetis vernus",
           "Centroptilum luteolum"	,
           "Chaetopteryx villosa"	,
           "Drusus annulatus"	,
           "Ephemerella ignita"	,
           "Ephemerella mucronata"	,
           "Isoperla goertzi"	,
           "Leuctra digitata"	,
           "Leuctra nigra"	,
           "Leuctra prima"	,
           "Nemoura cambrica"	,
           "Nemoura flexuosa"	,
           "Nemoura marginata"	,
           "Nemurella  pictetii"	,
           "Paraleptophlebia submarginata"	,
           "Plectrocnemia conspersa"	,
           "Potamophylax luctuosus"	,
           "Protonemura auberti"	,
           "Protonemura intricata"	,
           "Protonemura meyeri"	,
           "Rhyacophila fasciata"	,
           "Sericostoma personatum"	,
           "Silo pallipes"	,
           "Siphonoperla torrentium"	,
           "Tinodes rostocki"	,
           "Wormaldia occipitalis")

order1=c("T","P","T","P","E","E","E","T","T","E","E","P","P","P","P","P","P","P","P","E","T","T","P","P","P","T","T","T","P","T","T")
dft1$orders=orders
dft1 <- dft1[with(dft1, order(orders)),]
dft1=dft1[! dft1$phenophase %in% c("last", "first","value"), ]
taxa_com=cbind(taxaEPT,taxa1EPT)
taxa_com=data.frame(taxa_com)
taxa_com$spp=taxa_com$taxaEPT
taxa_com$order1=order1
dft1<-merge(taxa_com,dft1,by = c("spp"), all.x = TRUE, all.y = TRUE)

orders=c("T","T","P","P","T","T","E","E","E","E","P","P","E","E","T","T","T","T","E","E","E","E","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","E","E","T","T","T","T","P","P","P","P","P","P","T","T","T","T","T","T","P","P","T","T","T","T")

ord1 <- orders[seq(1, length(orders), 2)]
dft1$orders=orders


#dft1$taxa1EPT<- factor(dft1$taxa1EPT, levels = unique(dft1$taxa1EPT))
dft1 <- dft1[with(dft1, order(orders)),]
#dft1 <- dft1[with(dft1, order(-(as.factor(taxa1EPT))),]
#dft1x=dft1[order(dft1$taxa1EPT), ]
dft1<-with(dft1, dft1[order(orders,taxa1EPT),])
g1=ggplot(dft1, aes(taxa1EPT, day1)) + geom_point(size=10,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))+scale_size_manual(values =c("significant" = 6,"insignificant" = 4))+coord_flip()+ geom_hline(yintercept=0,lty=2)+ ylab("[shift, days per year]") +theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))            
g1t=g1+facet_wrap(~phenophase,scales = "free_x", ncol = 3, switch = "y") +ylab(NULL) +theme(strip.background = element_blank())+ geom_point(size=5,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))

g1t=g1t+theme_bw()
g1t=g1t+theme(legend.position="top")
g1t=g1t+theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))+ theme(text = element_text(size = 40))
italic.10.text <- element_text(face = "italic", size = 30)
g1t=g1t+theme(axis.text.y = italic.10.text)+coord_flip()+theme(strip.background = element_blank(), strip.text = element_blank())



#######################

#abundance change plot only
#######################################
#######################################
#######################################
taxa1EPT=c("Agapetus fuscipies"	,
           "Amphinemura standfussi"	,
           "Apatania fimbriata"	,
           "Brachyptera risi",
           "Baetis rohdani"	,
           "Baetis vernus",
           "Centroptilum luteolum"	,
           "Chaetopteryx villosa"	,
           "Drusus annulatus"	,
           "Ephemerella ignita"	,
           "Ephemerella mucronata"	,
           "Isoperla goertzi"	,
           "Leuctra digitata"	,
           "Leuctra nigra"	,
           "Leuctra prima"	,
           "Nemoura cambrica"	,
           "Nemoura flexuosa"	,
           "Nemoura marginata"	,
           "Nemurella  pictetii"	,
           "Paraleptophlebia submarginata"	,
           "Plectrocnemia conspersa"	,
           "Potamophylax luctuosus"	,
           "Protonemura auberti"	,
           "Protonemura intricata"	,
           "Protonemura meyeri"	,
           "Rhyacophila fasciata"	,
           "Sericostoma personatum"	,
           "Silo pallipes"	,
           "Siphonoperla torrentium"	,
           "Tinodes rostocki"	,
           "Wormaldia occipitalis")




a1$spp=a1$variable
dft2$phenophase = factor(dft2$phenophase, levels = c("value"))

#reordering columns to achive desirable position of the 
taxa_com=cbind(taxaEPT,taxa1EPT)
taxa_com=data.frame(taxa_com)
taxa_com$spp=taxa_com$taxaEPT
#dft2<-merge(taxa_com,dft2,by = c("spp"), all.x = TRUE, all.y = TRUE)
#dft2<-merge(taxa_com,dft2,by = c("spp"), all.x = TRUE, all.y = TRUE)
orders=c("T","T","P","P","T","T","E","E","E","E","P","P","E","E","T","T","T","T","E","E","E","E","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","E","E","T","T","T","T","P","P","P","P","P","P","T","T","T","T","T","T","P","P","T","T","T","T")

ord1 <- orders[seq(1, length(orders), 2)]
dft2$orders=ord1
#dft2=dft1[! dft1$phenophase %in% c("last", "first","ct","duation"), ]
dft2 <- dft2[with(dft2, order(orders)),]


#dft2$taxa1EPT <- factor(dft2$taxa1EPT, levels = unique(dft2$taxa1EPT))
#dft2 <- dft2[with(dft2, order(-taxa1EPT)),]

g2=ggplot(dft2, aes(taxa1EPT, day1)) + geom_point(size=10,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))+scale_size_manual(values =c("significant" = 6,"insignificant" = 4))+coord_flip()+ geom_hline(yintercept=0,lty=2)+ ylab("Changes in abundance [shift, days per year]") +theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))            
#g2t=g2+facet_wrap(~phenophase,scales = "free_x", ncol = 3, switch = "y", labeller = as_labeller(c(value="shift, abundance") ) ) +ylab(NULL) +theme(strip.background = element_blank())+ geom_point(size=5,aes(colour = as.factor(final)))+scale_color_manual(name = "significance",values = c("significant" = "#3399FF","insignificant" = "#333333"))

g2t=g2+theme_bw()
g2t=g2t+theme(legend.position="top")
g2t=g2t+theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))+ theme(text = element_text(size = 40))
italic.10.text <- element_text(face = "italic", size = 30)
g2t=g2t+ theme(axis.text.y = italic.10.text)+theme(strip.background = element_blank(), strip.text = element_blank())
#supplementary figure S3 

#total abundance code


#community level phenology
bn1=s_pheno[c(1,28:30,5,8,24:26)]
bn1=na.omit(b1)

#bn1$duration=bn1$week.last-bn1$week.first
bnx1=aggregate(bn1[, c(7:10)], list(bn1$Jahr,bn1$variable),mean,na.rm=TRUE, na.action=NULL) 
bnx2=aggregate(bnx1[, c(3:6)], list(bnx1$Group.1),mean,na.rm=TRUE, na.action=NULL) 
bnx3=aggregate(bnx1[, c(3:6)], list(bnx1$Group.1),mean,na.rm=TRUE, na.action=NULL) 

bnx=aggregate(bn1[, c(26:27)], list(bn1$Jahr),mean,na.rm=TRUE, na.action=NULL) 
bnx$year=bnx$Group.1
bnx=subset(bnx,bnx$year!="2006")
#mblm slopes of the phenophases
durat1=mblm(duration~year,data=bnx)
ct1=mblm(ct.week~year,data=bnx)


bnxm <- melt(bnx, id=c("Group.1"))
bn1$year=bn1$Group.1
bnx$year=bnx$Group.1
div.ab3<-merge(bnx,div.ab2, by = "year", all.x = TRUE, all.y = TRUE)
div.ab3$s_dur=scale(div.ab3$duration)+10
div.ab3$S_ct=scale(div.ab3$ct.week)+10

#model community level duration of emergence

#Duration <- gls(s_dur~mean_s+pattern+mean_s*pattern+sab,data=div.ab3, correlation=corAR1(), method="ML",na.action = na.omit)
#add diversity metrics to the model
Duration1 <- gls(s_dur~mean_s+pattern+mean_s*pattern+sab+S_rich+S_div,data=div.ab3, correlation=corAR1(), method="ML",na.action = na.omit)
#Duration2 <- gls(s_dur~(mean_s+pattern+sab+S_rich+S_div+turnover)*mean_s,data=div.ab3, correlation=corAR1(), method="ML",na.action = na.omit)
#plot duration
ggplot(div.ab3,aes(x=mean,y=duration,colour=as.factor(pattern)))+geom_point()+geom_smooth(method="lm")
#model community level Central tendency of emergence
#Peak_emergence <- gls(S_ct~(mean_s+pattern+mean_s*pattern+sab),data=div.ab3, correlation=corAR1(), method="ML",na.action = na.omit)
#plot CT 
ggplot(div.ab3,aes(x=mean,y=ct.week,colour=as.factor(pattern)))+geom_point()+geom_smooth(method="lm")

#model peak of emergencel
Peak_emergence1 <- gls(S_ct~(mean_s+pattern+mean_s*pattern+sab+S_rich+S_div),data=div.ab3, correlation=corAR1(), method="ML",na.action = na.omit)

cols <- c("week.first" = "#3366CC", "week.last" = "#33CC99", "ct.week" = "#CC99FF", "duration" = "#666666")
bnxm <- melt(bnx, id=c("Group.1","year"))
bnxm=bnxm[! bnxm$variable %in% c("week.first","week.last"), ]

p5x= ggplot(bnxm, aes(Group.1,value,shape=as.factor(variable)))+geom_point(size=5)+geom_smooth(method="lm",se = FALSE,lwd=1)+theme_bw()+ylab("Number of week")+xlab("Year")+theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))+geom_line(lty=2)
#shift in phenophases over time
p5x=p5x+ theme(legend.text = element_text(colour="black", size = 30, face = "bold"))
p5x=p5x+ theme(legend.position="top") #figure 4
#p5x=p5x
#p5x=p5x
temp=div.ab1[,c(1,24)]
temp=aggregate(temp[,2], list(temp$year),mean)#sum of all selected for phenology taxa
bnxm$year=bnxm$Group.1
temp$year=temp$Group.1####
bnxm<-merge(bnxm,temp, by = "Group.1", all.x = TRUE, all.y = TRUE)
bnx1<-merge(bnx,temp, by = "Group.1", all.x = TRUE, all.y = TRUE)
p5x=p5x+theme(legend.position = "none")
p5x=p5x+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
p5x=p5x+ theme(legend.position="none")
p5x=p5x+scale_color_manual(values = cols1)+ theme(legend.position="none")

