#Data on general abundance and diversity from trap B 

# subsetting and preparing Breitenbach datasets for the analysis#
#package  installation
#install.packages(c("plyr", "dplyr", "reshape2", "ggplot2", "nlme","gridExtra",
#"FD","mblm","trend","broom","ggfortify","graphics","Kendall","DataCombine",
#"modifiedmk","boot","vegan","gtools","codyn"))
require(plyr)
require(dplyr)
require(reshape2)
require(nlme)
require(car)
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
require(codyn)
# load package
library(sjPlot)
library(sjmisc)
library(sjlabelled)

###########################################
##########################################
#this is just my directory, plug yours here

#upload total abundance/ diversity matrix for trap B yrs 1969-2006
ept1= read.csv("ept_69_2006.csv",sep=";")
ept1[is.na(ept1)] <- 0 
taxa=ept1[,4:245]#separate taxa from the date columns

t1= read.table("bver1.txt",sep="\t",header = TRUE) #read Baetis vernus separately, since data are missing for reasons???
#from the main file, as column in the original file seems to be damaged, had to retype it
t1$II_Bae_ver_m1=as.numeric(t1$II_Bae_ver_m)
taxa$II_Bae_ver_m1=t1$II_Bae_ver_m1 #blend B.vernalis toogether with the rest
taxa=taxa[,2:243]
taxa=taxa[,c(242,1:241)]

#summ male and female abundance for every species
output <-data.frame( Bae_ver=apply(taxa[1:2], 1, sum),Bae_roh=apply(taxa[3:4], 1, sum), Bae_sp=apply(taxa[5:6], 1, sum),Eph_muc=apply(taxa[7:8],1, sum),
                     Eph_ign=apply(taxa[9:10],1, sum),Cen_lut=apply(taxa[11:12],1, sum),Par_sub=apply(taxa[13:14],1, sum),Ecd_sp=apply(taxa[15:16],1, sum),
                     Hab_fus=apply(taxa[17:18],1, sum),Hab_lau=apply(taxa[19:20],1, sum),Epe_sp=apply(taxa[21:22],1, sum), Bae_mut=apply(taxa[23:24],1, sum),
                     Bae_fus=apply(taxa[25:26],1, sum),Bae_sca=apply(taxa[27:28],1, sum),Eph_dan=apply(taxa[29:30],1, sum),Rhi_sp=apply(taxa[31:32],1, sum),
                     EEE=apply(taxa[33:35],1, sum),Amp_sta=apply(taxa[36:37],1, sum),Amp_sul=apply(taxa[38:39],1, sum),Bra_ris=apply(taxa[40:41],1, sum),
                     Bra_set=apply(taxa[42:43],1, sum),Iso_gra=apply(taxa[44:45],1, sum),Iso_goe=apply(taxa[46:47],1, sum),Leu_dig=apply(taxa[48:49],1, sum),
                     Leu_fus=apply(taxa[50:51],1, sum),Leu_ine=apply(taxa[52:53],1, sum),Leu_nig=apply(taxa[54:55],1, sum),Leu_pri=apply(taxa[56:57],1, sum),
                     Leu_sp=apply(taxa[58:59],1, sum),Nem_cam=apply(taxa[60:61],1, sum),Nem_cin=apply(taxa[62:63],1, sum),Nem_fle=apply(taxa[64:65],1, sum),
                     Nem_mar=apply(taxa[66:67],1, sum),Nrl_pic=apply(taxa[68:69],1, sum),Nem_sp=apply(taxa[70],1, sum),
                     Per_dis=apply(taxa[71:72],1, sum),Per_mic=apply(taxa[73:74],1, sum),Pro_aub=apply(taxa[75:76],1, sum),Pro_int=apply(taxa[77:78],1, sum),
                     Pro_mey=apply(taxa[79:80],1, sum),Pro_nit=apply(taxa[81:82],1, sum),Sip_tor=apply(taxa[83:84],1, sum),PPP=apply(taxa[85:87],1, sum)
                     ,Adi_fil=apply(taxa[88:89],1, sum),Adi_red=apply(taxa[90:91],1, sum),Aga_fus=apply(taxa[92:93],1, sum),Agr_mul=apply(taxa[94:95],1, sum)
                     ,Ana_ner=apply(taxa[96:97],1, sum),Ann_obs=apply(taxa[98:99],1, sum),Apa_fim=apply(taxa[100:101],1, sum),Ath_ate=apply(taxa[102:103],1, sum),
                     Ber_eid=apply(taxa[104:105],1, sum),Ber_pul=apply(taxa[106:107],1, sum), Bra_sub=apply(taxa[108:109],1, sum),Cer_dis=apply(taxa[110:111],1, sum),
                     Cha_vil=apply(taxa[112:113],1, sum),Che_lep=apply(taxa[114:115],1, sum),Cru_irr=apply(taxa[116:117],1, sum),Cyr_tri=apply(taxa[118:119],1, sum),
                     Dip_fel=apply(taxa[120:121],1, sum),Dru_ann=apply(taxa[122:123],1, sum),Dru_big=apply(taxa[124:125],1, sum),Glo_con=apply(taxa[126:127],1, sum),
                     Gly_pel=apply(taxa[128:129],1, sum), Gra_nig=apply(taxa[130:131],1, sum),Gra_nit=apply(taxa[132:133],1, sum),Hal_dig=apply(taxa[134:135],1, sum),
                     Hal_rad=apply(taxa[136:137],1, sum),Hal_tes=apply(taxa[138:139],1, sum),Hyd_ang=apply(taxa[140:141],1, sum),Hyd_ins=apply(taxa[142:143],1, sum),
                     Hyd_pel=apply(taxa[144:145],1, sum),Hyd_sax=apply(taxa[146:147],1, sum),Hyd_sil=apply(taxa[148:149],1, sum),Hyd_sp=apply(taxa[150:151],1, sum),
                     Hpt_ang=apply(taxa[152:153],1, sum),Hyd_pul=apply(taxa[154:155],1, sum),Hpt_spa=apply(taxa[156:157],1, sum),Las_bas=apply(taxa[158:159],1, sum),
                     Lim_aur=apply(taxa[160:161],1, sum),Lim_bip=apply(taxa[162:163],1, sum),Lim_cen=apply(taxa[164:165],1, sum),Lim_ext=apply(taxa[166:167],1, sum),
                     Lim_fus=apply(taxa[168:169],1, sum),Lim_gri=apply(taxa[170:171],1, sum),Lim_hir=apply(taxa[172:173],1, sum),Lim_lun=apply(taxa[174:175],1, sum),
                     Lim_ign=apply(taxa[176:177],1, sum),Lim_rho=apply(taxa[178:179],1, sum),Lim_spa=apply(taxa[180:181],1, sum),Lyp_pha=apply(taxa[182:183],1, sum),
                     Lip_red=apply(taxa[184:185],1, sum),Mic_lon=apply(taxa[186:187],1, sum),Mpt_lat=apply(taxa[188:189],1, sum),Mpt_nyc=apply(taxa[190:191],1, sum),
                     Mpt_seq=apply(taxa[192:193],1, sum),Mys_lon=apply(taxa[194:195],1, sum),Mys_nig=apply(taxa[196:197],1, sum),Odo_alb=apply(taxa[198:199],1, sum),
                     Par_pic=apply(taxa[200:201],1, sum),Ple_con=apply(taxa[202:203],1, sum),Pol_fla=apply(taxa[204:205],1, sum),Pot_cin=apply(taxa[206:207],1, sum),
                     Pot_lat=apply(taxa[208:209],1, sum),Pot_luc=apply(taxa[210:211],1, sum),Pot_nig=apply(taxa[212:213],1, sum),Psy_pus=apply(taxa[214:215],1, sum),
                     Pti_gra=apply(taxa[216:217],1, sum),Rhy_fas=apply(taxa[218:219],1, sum),Rhy_nub=apply(taxa[220:221],1, sum),Ser_fla=apply(taxa[222:223],1, sum),
                     Ser_per=apply(taxa[224:225],1, sum),Sil_pal=apply(taxa[226:227],1, sum),Ste_per=apply(taxa[228:229],1, sum),Syn_mos=apply(taxa[230:231],1, sum),
                     Tin_pal=apply(taxa[232:233],1, sum),Tin_ros=apply(taxa[234:235],1, sum),Tin_wae=apply(taxa[236:237],1, sum),Wor_occ=apply(taxa[238:239],1, sum),
                     TTT=apply(taxa[240:242],1, sum))
#separate date data
date=ept1[,1:3]
#merge taxa/ abundance and data
taxa1=cbind(date,taxa)#with species not summarized by sex
taxax=cbind(date,output)#with species summarized by sex

#melt dataset
eptmelt <- melt(taxa1, id=c("II_Jahr","II_Monat","II_Tag"))
#aggregate datsets by different grouping variables
eptx1=aggregate(eptmelt[, 5], list(eptmelt$II_Jahr,eptmelt$ variable),sum) #by year+ taxa

vec=c(rep("Eph",34),rep("Ple",52),rep("Tri",156))
#order names for every taxa in dataset
a1=data.frame(a1)
a2=cbind(a1,vec)
a2$Group.2=a2$a1
#eptx1<-merge(eptx1,a2, by = "Group.2", all.x = TRUE, all.y = TRUE)

#read in data for 2007-2010
ept2= read.csv("eptx2.csv",sep=";")

ept2[is.na(ept2)] <- 0 #replace NA with "0"=

taxa=ept2[,4:248]
output1<-data.frame( Bae_ver=apply(taxa[1:2], 1, sum),Bae_roh=apply(taxa[3:4], 1, sum), Bae_sp=apply(taxa[5:6], 1, sum),Eph_muc=apply(taxa[7:8],1, sum),
                     Eph_ign=apply(taxa[9:10],1, sum),Cen_lut=apply(taxa[11:12],1, sum),Par_sub=apply(taxa[13:14],1, sum),Ecd_sp=apply(taxa[15:16],1, sum),
                     Hab_fus=apply(taxa[17:18],1, sum),Hab_lau=apply(taxa[19:20],1, sum),Epe_sp=apply(taxa[21:22],1, sum), Bae_mut=apply(taxa[23:24],1, sum),
                     Bae_fus=apply(taxa[25:26],1, sum),Bae_sca=apply(taxa[27:28],1, sum),Eph_dan=apply(taxa[29:30],1, sum),Rhi_sp=apply(taxa[31:32],1, sum),
                     EEE=apply(taxa[33:35],1, sum),Amp_sta=apply(taxa[36:37],1, sum),Amp_sul=apply(taxa[38:39],1, sum),Bra_ris=apply(taxa[40:41],1, sum),
                     Bra_set=apply(taxa[42:43],1, sum),Iso_gra=apply(taxa[44:45],1, sum),Iso_goe=apply(taxa[46:47],1, sum),Leu_dig=apply(taxa[48:49],1, sum),
                     Leu_fus=apply(taxa[50:51],1, sum),Leu_ine=apply(taxa[52:53],1, sum),Leu_nig=apply(taxa[54:55],1, sum),Leu_pri=apply(taxa[56:57],1, sum),
                     Leu_sp=apply(taxa[58:60],1, sum),Nem_cam=apply(taxa[61:62],1, sum),Nem_cin=apply(taxa[63:64],1, sum),Nem_fle=apply(taxa[65:66],1, sum),
                     Nem_mar=apply(taxa[67:68],1, sum),Nrl_pic=apply(taxa[69:70],1, sum),Nem_sp=apply(taxa[71:72],1, sum),
                     Per_dis=apply(taxa[73:74],1, sum),Per_mic=apply(taxa[75:76],1, sum),Pro_aub=apply(taxa[77:78],1, sum),Pro_int=apply(taxa[79:80],1, sum),
                     Pro_mey=apply(taxa[81:82],1, sum),Pro_nit=apply(taxa[83:84],1, sum),Sip_tor=apply(taxa[85:86],1, sum),PPP=apply(taxa[87:90],1, sum)
                     ,Adi_fil=apply(taxa[91:92],1, sum),Adi_red=apply(taxa[93:94],1, sum),Aga_fus=apply(taxa[95:96],1, sum),Agr_mul=apply(taxa[97:98],1, sum)
                     ,Ana_ner=apply(taxa[99:100],1, sum),Ann_obs=apply(taxa[101:102],1, sum),Apa_fim=apply(taxa[103:104],1, sum),Ath_ate=apply(taxa[105:106],1, sum),
                     Ber_eid=apply(taxa[107:108],1, sum),Ber_pul=apply(taxa[109:110],1, sum), Bra_sub=apply(taxa[111:112],1, sum),Cer_dis=apply(taxa[113:114],1, sum),
                     Cha_vil=apply(taxa[115:116],1, sum),Che_lep=apply(taxa[117:118],1, sum),Cru_irr=apply(taxa[119:120],1, sum),Cyr_tri=apply(taxa[121:122],1, sum),
                     Dip_fel=apply(taxa[123:124],1, sum),Dru_ann=apply(taxa[125:126],1, sum),Dru_big=apply(taxa[127:128],1, sum),Glo_con=apply(taxa[129:130],1, sum),
                     Gly_pel=apply(taxa[131:132],1, sum), Gra_nig=apply(taxa[133:134],1, sum),Gra_nit=apply(taxa[135:136],1, sum),Hal_dig=apply(taxa[137:138],1, sum),
                     Hal_rad=apply(taxa[139:140],1, sum),Hal_tes=apply(taxa[141:142],1, sum),Hyd_ang=apply(taxa[143:144],1, sum),Hyd_ins=apply(taxa[145:146],1, sum),
                     Hyd_pel=apply(taxa[147:148],1, sum),Hyd_sax=apply(taxa[149:150],1, sum),Hyd_sil=apply(taxa[151:152],1, sum),Hyd_sp=apply(taxa[153:154],1, sum),
                     Hpt_ang=apply(taxa[155:156],1, sum),Hyd_pul=apply(taxa[157:158],1, sum),Hpt_spa=apply(taxa[159:160],1, sum),Las_bas=apply(taxa[161:162],1, sum),
                     Lim_aur=apply(taxa[163:164],1, sum),Lim_bip=apply(taxa[165:166],1, sum),Lim_cen=apply(taxa[167:168],1, sum),Lim_ext=apply(taxa[169:170],1, sum),
                     Lim_fus=apply(taxa[171:172],1, sum),Lim_gri=apply(taxa[173:174],1, sum),Lim_hir=apply(taxa[175:176],1, sum),Lim_lun=apply(taxa[177:178],1, sum),
                     Lim_ign=apply(taxa[179:180],1, sum),Lim_rho=apply(taxa[181:182],1, sum),Lim_spa=apply(taxa[183:184],1, sum),Lyp_pha=apply(taxa[185:186],1, sum),
                     Lip_red=apply(taxa[187:188],1, sum),Mic_lon=apply(taxa[189:190],1, sum),Mpt_lat=apply(taxa[191:192],1, sum),Mpt_nyc=apply(taxa[193:194],1, sum),
                     Mpt_seq=apply(taxa[195:196],1, sum),Mys_lon=apply(taxa[197:198],1, sum),Mys_nig=apply(taxa[199:200],1, sum),Odo_alb=apply(taxa[201:202],1, sum),
                     Par_pic=apply(taxa[203:204],1, sum),Ple_con=apply(taxa[205:206],1, sum),Pol_fla=apply(taxa[207:208],1, sum),Pot_cin=apply(taxa[209:210],1, sum),
                     Pot_lat=apply(taxa[211:212],1, sum),Pot_luc=apply(taxa[213:214],1, sum),Pot_nig=apply(taxa[215:216],1, sum),Psy_pus=apply(taxa[217:218],1, sum),
                     Pti_gra=apply(taxa[219:220],1, sum),Rhy_fas=apply(taxa[221:222],1, sum),Rhy_nub=apply(taxa[223:224],1, sum),Ser_fla=apply(taxa[225:226],1, sum),
                     Ser_per=apply(taxa[227:228],1, sum),Sil_pal=apply(taxa[229:230],1, sum),Ste_per=apply(taxa[231:232],1, sum),Syn_mos=apply(taxa[233:234],1, sum),
                     Tin_pal=apply(taxa[235:236],1, sum),Tin_ros=apply(taxa[237:238],1, sum),Tin_wae=apply(taxa[239:240],1, sum),Wor_occ=apply(taxa[241:242],1, sum),
                     TTT=apply(taxa[243:245],1, sum))
out=rbind(output,output1)
date1=ept2[,1:3]
taxay=cbind(date1,output1)#2007-2010

taxaf=rbind(taxax,taxay)#combined sex summarized taxa for 1969-2005 and 2007-2010
taxaf<- subset(taxaf, II_Jahr != c("0"))#remove empty year rows
taxaf<- subset(taxaf, II_Jahr != c("2006"))#remove 2006

ept2m <- melt(taxaf, id=c("II_Jahr","II_Monat","II_Tag"))#melt taxa set 1969-2005 and 2007-2010


sp=read.table("spp.txt",sep="\t",header = TRUE) #reading in taxa full names
#total taxa/abundance 1969- /2006/-2010
ml<-merge(sp,ept2m, by = "variable", all.x = TRUE, all.y = TRUE)
################################################
#################################################
###################################################
#Analysis of the trophic groups relative abundance

eptw <- dcast(ml,Taxon~II_Jahr, value.var="value", fun=sum)# spp level abundance matrix in wide form
#"eptw" - eptw- "ephemeroptera-plecoptera-trichoptera wide"
mtr= read.table("metr.txt",sep="\t",header = TRUE) #reading in metrics of trophic abundance
#calculated on the eptw using Asterics software
troph=mtr[,c(1,46,49:50,52:53)]#removing unused metrics
troph<- subset(troph, Year != c("2006"))#remove 2006
#logit transfor trophic groups relative abundance
troph1=troph
troph1[is.na(troph1)] <- 0

troph1$logit_gra=log(troph1$X......Grazers.and.scrapers)+10
troph1$logit_shre=log(troph1$X......Shredders)+10
troph1$logit_gat=log(troph1$X......Gatherers.Collectors)+10
troph1$logit_filt=log(troph1$X......Passive.filter.feeders)+10
troph1$logit_pre=log(troph1$X......Predators)+10
troph1[is.na(troph1)] <- 0
troph1$year=troph1$Year
#sen's slopes of the change in relative abunce
##of the trophic groups
summary(mblm(logit_gra~Year,data=troph1))
graz=mblm(logit_gra~year,data=troph1)
summary(mblm(logit_shre~Year,data=troph1))
shred=mblm(logit_shre~year,data=troph1)
summary(mblm(logit_gat~Year,data=troph1))
gath=mblm(logit_gat~year,data=troph1)
summary(mblm(logit_filt~Year,data=troph1))
filtr=mblm(logit_filt~year,data=troph1)
summary(mblm(logit_pre~Year,data=troph1))
pred=mblm(logit_pre~year,data=troph1)
#tab_model(graz,shred,gath,filtr,pred,p.style="a")

trp2m <- melt(troph, id=c("Year"))

########################################
#######################################
#######################################
#modelling community variables with environmental variables
#rad in dischrge data from Breitenbach
dis= read.table("discharge.txt",sep="\t",header = TRUE) #reading in abiotic [discharge] data
#rename
dis$year=dis$Jahr
#read temperature data
discharge1=dis[,c(1:13,15,16)]
#read_pattern of discharge
dp= read.csv("disch_pattern.csv",sep=",")#read in temperature data (water)
dp$mean_minus_sd=dp$mean-dp$sd #calculate mean discharge- sd
dp$mean_plus_sd=dp$mean+dp$sd #calculate mean discharge+ sd
dp2=dp[,c(1,2,4:6)]#delete sd as not used in the plot
#plot fig s2
dp1<- melt(dp2, id=c("var","pattern"))#melt to make ggplot possible
dp1$var <- factor(dp1$var,levels = c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep"))
plot_s2=ggplot(dp1,aes(x=as.factor(var),y=value,colour=as.factor(variable), group=as.factor(variable)))+geom_point()+geom_line(aes(colour=as.factor(variable)))+facet_wrap(~pattern)

plot_s2=plot_s2+theme_bw()+xlab("months")+ylab("monthly mean discharge, L*s-1")

#######################################################################
temp= read.csv("WaterTBreitenbach_final.csv",sep=",")#read in temperature data (water)
temperature1=mblm(AnnualMean~year,data=temp)
#merge discharge and temperature data
vars<-merge(temp,dis, by = "year", all.x = TRUE, all.y = TRUE)
air= read.csv("Breitenbach_AirTemp.csv",sep=",")#air temperature data
air=air[20:61,]#extract required yearsplot(temp$AnnualMean~temp$year,pch=19,ylim=c(6,10))
#points(air$AirT~air$Year)
#lines(temp$AnnualMean~temp$year,cex=1)
#lines(air$AirT~air$Year,cex=3,col="blue")

ins1=aggregate(ept2m[,5], list(ept2m$II_Jahr),sum)#general abundance of all insects
ins1$year=ins1$Group.1#rename column
write.table(ins1, "totab.txt", sep=",")#total abundance siteB


ins2<-merge(ins1,vars, by = "year", all.x = TRUE, all.y = TRUE)
#merge diversity and abundance data with environmental data
#rename for convinience
data=ins2 #rename
data<- subset(data, year != c("2006"))
#standartize temperature and abundance data
data$mean_s=scale(data$AnnualMean)+10
data$sab=scale(data$x)+10
data$abund1=(data$x)/1000
data$temperature=data$mean_s


attach(data)
#modelling abundance
abundance<-gls(sab~mean_s+pattern+mean_s*pattern,data=data,na.action=na.omit)
#summary table



#modelling factors influencing rlative abundance of the trophic groups
trp2m$year=trp2m$Year#rename column
vars1<-merge(data,trp2m, by = "year", all.x = TRUE, all.y = TRUE)#merge abundance data, enviro data with trophic metrics
#vars1=vars1[2:206,]
vars1$mean_s=scale(vars1$AnnualMean)+10
vars1$lvalue=logit(vars1$value*0.01) #logit relative troph group abundance
vars1$sval=scale(vars1$lvalue)+10#standartize relative troph group abundance
trophic_groups<- gls(sval~(mean_s+pattern+mean_s*pattern)*variable,data=vars1,na.action=na.omit)#modell for troph groups

#ptrop_time=ggplot(vars1,aes(x=mean_s,y=value, colour=as.factor(pattern)))+geom_point()+geom_smooth(method="lm",se=FALSE)+facet_wrap(~variable,scales="free")

#ggplot(vars1,aes(x=year,y=value, colour=as.factor(variable)))+geom_point()+geom_smooth(se=FALSE)

#analysing diversity using vegan package
#df1=smartbind(taxa1, ept2)
#df1<- subset(df1, II_Jahr != c("2006"))
#df1[is.na(df1)] <- 0
#tax=df1[,4:269]#separate taxa abundance matrix from the rest of the data
taxaf1=taxaf #to preserve dataset for the decade turnover calculation
taxaf=taxaf[,c(1,4:123)]
#taxaf=taxaf(2:43,)
taxaf=aggregate(taxaf[, 2:121], list(taxaf$II_Jahr),sum)
eptdiv=taxaf
tax=taxaf[,2:121]

H <- diversity(tax)#Shannon's diversity
simp <- diversity(tax, "simpson")#simpson diversity
invsimp <- diversity(tax, "inv")#simpsom invariable
S <- specnumber(tax) #species number/ richness
J <- H/log(S)#evennes
#bind diversity metrics
div_ept=cbind(H,S,J)# dataframe with diversity indices
div_ept=data.frame(div_ept)#as data frame
div_ept$habitat=rep()#????
div_ept$year=c(1969:2005,2007:2010)#remove 2006

#combine taxa matrix with diversity indeces
taxaf$richness=S #species richness
taxaf$div=H #shannons diversity
taxaf$simp1=simp #Simson's diversity
taxaf$evennes=J #evennes
#calculate sen's slope models
#summary(mblm(evennes~Group.1,data=taxaf))#evennes
div=taxaf[,c(1,122:125)]
ab=aggregate(ept2m[, 5], list(ept2m$II_Jahr),sum)
ab$II_Jahr=ab$Group.1
div$II_Jahr=div$Group.1
div.ab<-merge(div,ab, by = "II_Jahr", all.x = TRUE, all.y = TRUE)
div.ab$year=div.ab$II_Jahr
div.ab1<-merge(div.ab,vars1, by = "year", all.x = TRUE, all.y = TRUE)
#scale the diversity indeces
div.ab1$s_mean=scale(div.ab1$AnnualMean)+10
div.ab1$S_rich=scale(div.ab1$richness)+10
div.ab1$S_div=scale(div.ab1$div)+10
div.ab1$S_simp=scale(div.ab1$simp1)+10
div.ab1$S_ev=scale(div.ab1$evennes)+10
div.ab1$S_troph=scale(div.ab1$value)+10
attach(div.ab1)

#model diversity matricies
Shannos_diversity<- gls(S_div~mean_s+pattern+mean_s*pattern,data=div.ab1, correlation=corAR1(), method="ML",na.action = na.omit)

Species_richness<- gls(S_rich~mean_s+pattern+mean_s*pattern,data=div.ab1, correlation=corAR1(), method="ML",na.action = na.omit)

Species_evennes<- gls(S_ev~mean_s+pattern+mean_s*pattern,data=div.ab1, correlation=corAR1(), method="ML",na.action = na.omit)
#dis$Year=dis$Jahr
#trp<-merge(dis,troph, by = "Year", all.x = TRUE, all.y = TRUE)
ab1=aggregate(ept2m[, 5], list(ept2m$II_Jahr,ept2m$variable),sum)

#calculate community turnover using codyn package
total.res <- turnover(df=ab1,time.var = "Group.1",species.var = "Group.2",abundance.var = "x",replicate.var=NA)
rate=rate_change_interval(ab1,time.var = "Group.1", species.var = "Group.2", abundance.var = "x",replicate.var =NA)

#appearence component
total.ap <- turnover(df=ab1,time.var = "Group.1", species.var = "Group.2",abundance.var = "x",replicate.var=NA,metric="appearance")

#disappearence component
total.ds <- turnover(df=ab1,time.var = "Group.1", species.var = "Group.2",abundance.var = "x",replicate.var=NA,metric="disappearance")
total.res$cat=rep("total",40)

total.ap$cat=rep("appearence",40)
total.ds$cat=rep("disappearence",40)

total.res$year=total.res$Group.1
total.ap$year=total.ap$Group.1
total.ds$year=total.ds$Group.1


total.res$turnover=total.res$total
total.res=total.res[,2:5]

total.ap$turnover=total.ap$appearance
total.ap=total.ap[,2:5]

total.ds$turnover=total.ds$disappearance 
total.ds=total.ds[,2:5]

#total turnover data frame
tot.tr=rbind(total.res,total.ap)
tottr=rbind(tot.tr,total.ds)
tottr=subset(tottr,tottr$year!="2006")

#merge taxa abundance matrix with turnover df
div.ab2<-merge(div.ab1,total.res, by = c("year","Group.1"), all.x = TRUE, all.y = TRUE)
div.ab2=subset(div.ab2,div.ab2$year!="2006")
#########################################

#plotting things...

trp2m$year=trp2m$Year
tropm<-merge(trp2m,dis, by = "year", all.x = TRUE, all.y = TRUE)
#plotting dischrge data
pdis=ggplot(vars1,aes(x=year,y=mean.disch))+geom_point(aes(color = factor(pattern),size=5))+geom_line(lty=2, size=0.3) #discharge plot
pdis=pdis+theme_bw()+xlab("")+ylab(expression(bold(paste(~ 'Daily maximum discharge [L*s'^-1~ ']'))))
pdis=pdis+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
pdis=pdis+ theme(legend.position="none")
pdis
ggsave("Fig.1b_temp_time.png", width = 15, height = 15, units = "cm") 
#plot temperature
ptemp=ggplot(vars1,aes(x=year,y=AnnualMean))+geom_point(aes(size=5))+geom_smooth(method="lm",se=FALSE)+geom_line(lty=2, size=0.3)#temperature plot
ptemp=ptemp+theme_bw()+xlab("")+ylab("Water temperature [°C]")+ theme(legend.position="none")
ptemp=ptemp+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))
ptemp=ptemp+ theme(legend.position="none")
ptemp
ggsave("Fig.1a_temp_time.png", width = 15, height = 15, units = "cm") 
#plot abundance
pab=ggplot(data,aes(x=year,y=abund1))+geom_point(aes(size=5))+geom_smooth(method="lm",se=FALSE)+geom_line(lty=2, size=0.3) #abundance vs yrs plot
pab=pab+theme_bw()+xlab("")+ylab("Abundance [Individuals × 1000]")+ theme(legend.position="none")
pab=pab+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
pab=pab+ theme(legend.position="none")
pab
ggsave("Fig.2a_abundance_time.png", width = 15, height = 15, units = "cm") 
#plot abundance vs temperature
pabt=ggplot(data,aes(x=AnnualMean,y=x))+geom_point(aes(size=5))+geom_smooth(method="lm",se=FALSE)+geom_line(lty=2) #abundance vs temp plot
pabt=pabt+theme_bw()+xlab("mean annual temperature, °C")+ylab("summ of animal abundance, specimens")+ theme(legend.position="none")
pabt=pabt+theme(axis.text=element_text(size=14),axis.title=element_text(size=20,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))

grid.arrange(ptemp,pdis,ncol=2) #figure 1
#data=subset(data,data$year!="2006")
summary(mblm(x~year,data=data))# sens slope for abundance from year
data$abundance=data$x
summ=mblm(abundance~year,data=data)
#summary(mblm(x~AnnualMean,data=data1))# sens slope for abundance from t°C

#######################################
par(mfrow=c(2,2))
####################################################
#subset undersampled 2006 from diversity matrix
taxaf<- subset(taxaf, Group.1 != "2006")
#plot Shannons diversity over time
pdiv=ggplot(taxaf,aes(x=Group.1,y=div))+geom_point(aes(size=5))+geom_smooth(method="lm",se=FALSE)+geom_line(lty=2, size=0.3) #temperature plot
pdiv=pdiv+theme_bw()+xlab("Year")+ylab("Shannon diversity")+ theme(legend.position="none")
pdiv=pdiv+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
pdiv=pdiv+ theme(legend.position="none")

#plot Species richness over time
prich=ggplot(taxaf,aes(x=Group.1,y=richness))+geom_point(aes(size=5))+geom_smooth(method="lm",se=FALSE)+geom_line(lty=2, size=0.3) #temperature plot
prich=prich+theme_bw()+xlab("Year")+ylab("Species richness")+ theme(legend.position="none")
prich=prich+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
prich=prich+ theme(legend.position="none")

#plot evennes over time
pev=ggplot(taxaf,aes(x=Group.1,y=evennes))+geom_point(aes(size=5))+geom_smooth(method="lm",se=FALSE)+geom_line(lty=2, size=0.3) #temperature plot
pev=pev+theme_bw()+xlab("years")+ylab("Pielou evennes")+ theme(legend.position="none")
pev=pev+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
pev=pev+ theme(legend.position="none")
#arrange grid


#plot evennes
pev=ggplot(taxaf,aes(x=Group.1,y=evennes))+geom_point(aes(size=5))+geom_smooth(method="lm",se=FALSE)+geom_smooth(colour=("darkred"),se=FALSE)+geom_line(lty=2, size=0.3) #temperature plot
pev=pev+theme_bw()+xlab("Year")+ylab("Species evenness")+ theme(legend.position="none")
pev=pev+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
pev=pev+ theme(legend.position="none")
grid.arrange(pab,prich,pdiv,pev,ncol=2)#figure2
ggsave("Fig.2.png", width = 30, height = 30, units = "cm") 





taxaf[is.na(taxaf)] <- 0#
taxaf=subset(taxaf,taxaf$year!="2006")
#compute sen's slopes for the community diversity metricies
summary(mblm(S_rich~year,data=div.ab1))
rich=mblm(S_rich~year,data=div.ab1)
summary(mblm(S_div~year,data=div.ab1))
div=mblm(S_div~year,data=div.ab1)
summary(mblm(S_simp~Group.1,div.ab1))
summary(mblm(S_ev~year,data=div.ab1))
evenne=mblm(S_ev~year,data=div.ab1)

#plot trophic groups relative abundance 
cols <- c("X......Grazers.and.scrapers" = "#009999", "X......Shredders" = "#330099", "X......Gatherers.Collectors" = "#666666", "X......Passive.filter.feeders" = "#CC6600","X......Predators"="#CC3399")
ptroph=ggplot(div.ab2,aes(x=year,y=value, colour=as.factor(variable)))+geom_point(size=3)+geom_smooth(method="lm",se=FALSE)
ptroph=ptroph+theme_bw()+xlab("Year")+ylab("Relative abundance [%]")
ptroph=ptroph+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
ptroph=ptroph+ theme(legend.position="none")
ptroph=ptroph+scale_color_manual(values = cols)+ theme(legend.position="none")
ptroph
ggsave("Fig.3a.png", width = 15, height = 15, units = "cm") 

#####################################################
h1=data.frame(H)
h1$year=c(1969:2005,2007:2010)
div.ab2<-merge(h1,div.ab2, by = "year", all.x = TRUE, all.y = TRUE)

#plot turnover over time

cols1 <- c("appearence" = "#95c11f", "disappearence" = "darkred", "total" = "darkgrey")
turn=ggplot(tottr,aes(x=year,y=turnover,colour=as.factor(cat)))+geom_point(size=3)+geom_smooth(method="lm",se=FALSE)
turn=turn+theme_bw()+xlab("Year")+ylab("Turnover")
turn=turn+theme(axis.text=element_text(size=14),axis.title=element_text(size=20,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
turn=turn+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
turn=turn+ theme(legend.position="none")
turn=turn+scale_color_manual(values = cols1)+ theme(legend.position="none")
turn_p=turn

turn
ggsave("Fig.3b.png", width = 15, height = 15, units = "cm") 
#abundance through the all sites
#plot abundance over time for all available traps
colsx <- c("A" = "#666600", "B" = "#FF3300", "C" = "#663399")
abtot= read.table("all sites abund.txt",sep="\t",header = TRUE) #reading in river data
abtot<- subset(abtot,Group.1 != "2006")
pa=ggplot(abtot,aes(Group.1,x,colour=as.factor(label)))+geom_point()+geom_line()+theme_bw()
pa=pa+scale_color_manual(values = colsx)+ theme(legend.position="none")
pa+xlab("years")+ylab("abundance, specimens")
abtotA=abtot[abtot$label == "B" | abtot$label == "A", ]
abtotA=subset(abtotA, Group.1>1983 & Group.1<1996)
abtotC=abtot[abtot$label == "B" | abtot$label == "C", ]
abtotC=subset(abtotC, Group.1>1974 & Group.1<1986)

#plot traps A,B, C together
abtot1<- subset(abtot, label != "E")
abtot1<- subset(abtot1, label != "G")
pa1=ggplot(abtot1,aes(Group.1,x,colour=as.factor(label)))+geom_point(size=5)+geom_line(lwd=2)+theme_bw()
pa1=pa1+scale_color_manual(values = colsx)+ theme(legend.position="none")
pa1+xlab("years")+ylab("abundance, specimens")+theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))+ theme(text = element_text(size = 35))
pa1=pa1+scale_color_manual(values = colsx)+ theme(legend.position="none")
pa1=pa1+xlab("years")+ylab("abundance, specimens")
colsA <- c( "B" = "#FF3300", "A" = "#666600")
pA=ggplot(abtotA,aes(Group.1,x,colour=as.factor(label)))+geom_point()+geom_line()+theme_bw()
pA=pA+scale_color_manual(values = colsA)+ theme(legend.position="none")
pA=pA+xlab("years")+ylab("abundance, specimens")+xlim(1975, 1996)+ylim(0,40000)
colsC <- c( "B" = "#FF3300", "C" = "#663399")
pC=ggplot(abtotC,aes(Group.1,x,colour=as.factor(label)))+geom_point()+geom_line()+theme_bw()
pC=pC+scale_color_manual(values = colsC)+ theme(legend.position="none")
pC=pC+xlab("years")+ylab("abundance, specimens")+xlim(1975, 1996)+ylim(0,40000)
#grid.arrange(pC,pA,ncol=2)

sA<- dcast(abtotA,Group.1~label, value.var="x", fun=sum)
sC<- dcast(abtotC,Group.1~label, value.var="x", fun=sum)

cor.test(sA$A,sA$B,method="kendall")
cor.test(sC$C,sC$B,method="kendall")

summary(mblm(x~year,data=data))# sens slope for abundance from year

sumx=mblm(x~year,data=data)
#######################################
####################################################

taxaf<- subset(taxaf, year != "2006")
#publication quality figures of the abundance and diversity indices
pab=ggplot(data,aes(x=year,y=abund1))+geom_point(size=15)+geom_smooth(method="lm",se=FALSE)+geom_line() #abundance vs yrs plot
pab=pab+theme_bw()+xlab("year")+ylab("Abundance, thouthands of specimens")+ theme(legend.position="none")
pab=pab+theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))+theme(legend.text = element_text(colour="black", size = 40, face = "bold"))
pdiv=ggplot(taxaf,aes(x=Group.1,y=div))+geom_point(size=15)+geom_smooth(method="lm",se=FALSE)+geom_line() #temperature plot
pdiv=pdiv+theme_bw()+xlab("year")+ylab("Shannon's diversity")+ theme(legend.position="none")
pdiv=pdiv+theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))+theme(legend.text = element_text(colour="black", size = 40, face = "bold"))
prich=ggplot(taxaf,aes(x=Group.1,y=richness))+geom_point(size=15)+geom_smooth(method="lm",se=FALSE)+geom_line() #temperature plot
prich=prich+theme_bw()+xlab("year")+ylab("Species richness")+ theme(legend.position="none")
prich=prich+theme(axis.text=element_text(size=40),axis.title=element_text(size=40,face="bold"))+theme(legend.text = element_text(colour="black", size = 40, face = "bold"))

#total abundance code
div.ab2=subset(div.ab2,div.ab2$year!="2006")#
#publication quality figures
cols <- c("X......Grazers.and.scrapers" = "#009999", "X......Shredders" = "#330099", "X......Gatherers.Collectors" = "#666666", "X......Passive.filter.feeders" = "#CC6600","X......Predators"="#CC3399")
ptroph=ggplot(div.ab2,aes(x=year,y=value, colour=as.factor(variable)))+geom_point(size=16)+geom_smooth(method="lm",se=FALSE)
ptroph=ptroph+theme_bw()+xlab("year")+ylab("Trophic groups share, %")
ptroph=ptroph+theme(axis.text=element_text(size=14),axis.title=element_text(size=20,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
pab=pab+theme(legend.text=element_text(size=15))+ theme(legend.position="none")

ptroph=ptroph+scale_color_manual(values = cols)+ theme(legend.position="none")

cols1 <- c("total" = "#9933FF")
tot1=subset(tottr,tottr$cat=="total")
#arrange figure display 
turn1=ggplot(tot1,aes(x=year,y=turnover,colour=as.factor(cat)))+geom_point(size=5)+geom_smooth(method="lm",se=FALSE)
turn1=turn1+theme_bw()+xlab("Year")+ylab("Turnover")
turn1=turn1+theme(axis.text=element_text(size=14),axis.title=element_text(size=20,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
turn1=turn1+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
turn1=turn1+ theme(legend.position="none")
turn1=turn1+scale_color_manual(values = cols1)+ theme(legend.position="none")
turnover. <- gls(turnover~mean_s+pattern+mean_s*pattern+sab+S_rich+S_div,data=div.ab2, correlation=corAR1(), method="ML",na.action = na.omit)

#interdecadal turnover
#check extinction and immigration
tax=taxaf1[,c(1,4:123)]
eptx1<- melt(tax, id=c("II_Jahr"))

#aggregate taxa abundance for each year
eptx1=aggregate(eptx1[, 3], list(eptx1$II_Jahr,eptx1$variable),sum)
eptx1$II_Jahr=eptx1$Group.1
dec1=subset(eptx1, eptx1$II_Jahr>1969 & eptx1$II_Jahr<1980)
dec2=subset(eptx1, eptx1$II_Jahr>1979 & eptx1$II_Jahr<1990)
dec3=subset(eptx1, eptx1$II_Jahr>1989 & eptx1$II_Jahr<2000)
dec4=subset(eptx1, eptx1$II_Jahr>1999 & eptx1$II_Jahr<2011)
dec1$decade=rep("first",1200)
dec2$decade=rep("second",1200)
dec3$decade=rep("third",1200)
dec4$decade=rep("last",1200)
decs=rbind(dec1,dec2,dec3,dec4)
#decs=aggregate(decs[, 3], list(decs$variable,decs$de),sum)



write.table(decs, file = "spp_decades.csv")
unique(eptx_1$Group.2)
interval<-c("1_2","1_2","2_3","2_3","3_4","3_4","first-last","first-last")
type=c("main","main","main","main","main","main","total","total")
direc=rep(c("gains","losses"),4)
spp=c(16,6,13,10,11,10,20,6)
turnover=cbind(interval,direc,spp,type)
turn=data.frame(turnover)
x=as.numeric(spp)
turn$spp=x
#ttot=subset(turn,turn$type=="total")
#tmain=subset(turn,turn$type=="main")
#hist1=ggplot(tmain,aes(x=interval,y=spp,fill=as.factor(direc)))+geom_bar(stat="identity",position=position_dodge())+
  #theme_minimal()+scale_fill_manual(values=c("#95c11f", "darkred"))

hist1=ggplot(turn,aes(x=interval,y=spp,fill=as.factor(direc)))+geom_bar(stat="identity",position=position_dodge())+
theme_minimal()+scale_fill_manual(values=c("#95c11f", "darkred"))
hist1=hist1+theme_bw()+xlab("Compared decade")+ylab("Number of species")
hist1=hist1+theme(axis.text=element_text(size=14),axis.title=element_text(size=20,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
hist1=hist1+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
hist1=hist1+theme(legend.position="none")
hist1=hist1+scale_color_manual(values = cols1)+ theme(legend.position="none")

hist1
ggsave("Fig.3c.png", width = 15, height = 15, units = "cm") 

hist2=ggplot(ttot,aes(x=interval,y=spp,fill=as.factor(direc)))+geom_bar(stat="identity",position=position_dodge())+
  theme_minimal()+scale_fill_manual(values=c("#95c11f", "darkred"))
hist2=hist2+theme_bw()+xlab("Compared decade")+ylab("Number of species")
hist2=hist2+theme(axis.text=element_text(size=14),axis.title=element_text(size=20,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
hist2=hist2+theme(axis.text=element_text(size=20),axis.title=element_text(size=18,face="bold"))+theme(legend.text = element_text(colour="black", size = 20, face = "bold"))#+ theme(legend.position="none")
hist2=hist2+theme(legend.position="none")
hist2=hist2+scale_color_manual(values = cols1)+ theme(legend.position="none")

hist2
ggsave("Fig.3d.png", width = 15, height = 15, units = "cm") 

grid.arrange(ptroph,turn_p,hist1,hist2,ncol=2)#figure 3
#table for model's output

tab_model(abundance, trophic_groups,Shannos_diversity,Species_richness,Species_evennes,turnover.,Duration1,Peak_emergence1,p.style="a")
tab_model(temperature1,summ,div,rich,evenne,ct1,durat1,graz,shred,gath,filtr,pred,p.style="a")
