
# 1) Download file "tg_0.25deg_reg_v17.0.nc.gz" from: https://www.ecad.eu
# 2) Extract .gz compressed files using:
#    library(R.utils)
#    R.utils::gunzip("tg_0.25deg_reg_v17.0.nc.gz")

# 3) Extract data:

library(sp)
library(raster)
library(ncdf4)

#Set coordinates of target point: 50°39'N 9°38'E
lon <- 9.633333  # longitude of location 					
lat <- 50.650000# latitude  of location 			

#Extract data, see:
#http://geog.uoregon.edu/GeogR/topics/netCDF-read-ncdf4.html
#https://stackoverflow.com/questions/20621200/extract-time-series-of-a-point-lon-lat-from-netcdf-in-r

#Mean daily temperature
ncin <- nc_open("tg_0.25deg_reg_v17.0.nc" )
print(ncin)
t <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(t)
nt
obsoutput <- ncvar_get(ncin, 
                       start= c(which.min(abs(ncin$dim$longitude$vals -   lon)), # look for closest long
                                which.min(abs(ncin$dim$latitude$vals -  lat)),  # look for closest lat
                                1),
                       count=c(1,1,-1))
DataMeanT <- data.frame(DateN= t, MeanDailyT = obsoutput)
nc_close(ncin)
head(DataMeanT)
#check if there are NAs =999
summary(DataMeanT)


Data=DataMeanT
Data$Date=as.Date(Data$DateN,origin="1950-01-01")
Data$Year=format(Data$Date,"%Y")
Data$Month=format(Data$Date,"%m")
head(Data)
Data$YearMonth=format(Data$Date, format="%Y-%b")

Data_annual=aggregate(("T_AnnualMean" = MeanDailyT) ~ Year,data = Data, FUN = mean,na.action = na.pass)
names(Data_annual)[2]<-"AirT"
head(Data_annual)

#Export table
write.table(Data_annual, "Breitenbach_AirTemp.csv", row.names = FALSE, append = FALSE, col.names = TRUE, sep = ", ", quote = TRUE)
