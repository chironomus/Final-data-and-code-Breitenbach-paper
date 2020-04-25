Temp= read.csv("Initial_Table_Temp.csv",sep=";")
head(Temp)
names(Temp)

# create row for year 2006  -----------------------------------------------
new.row <- head(Temp[NA,], 1)
new.row["year"] <- 2006
Temp=rbind(Temp,new.row)
Temp=Temp[order(Temp$year),]
row.names(Temp)=Temp$year

# Replace missing data with average of data from previous and next year -------------
# Except for year 1969 (first one of the time series) for which next the 2 years were used. 

Temp["1969","Jan"]<- mean(c(Temp["1970","Jan"],Temp["1971","Jan"]))
Temp["1969","Feb"]<- mean(c(Temp["1970","Feb"],Temp["1971","Feb"]))
Temp["1969","Mar"]<- mean(c(Temp["1970","Mar"],Temp["1971","Mar"]))
Temp["1984","Sep"]<- mean(c(Temp["1983","Sep"],Temp["1985","Sep"]))
Temp["1984","Oct"]<- mean(c(Temp["1983","Oct"],Temp["1985","Oct"]))
Temp["1984","Nov"]<- mean(c(Temp["1983","Nov"],Temp["1985","Nov"]))
Temp["1984","Dec"]<- mean(c(Temp["1983","Dec"],Temp["1985","Dec"]))
Temp["1985","Jan"]<- mean(c(Temp["1984","Jan"],Temp["1986","Jan"]))
Temp["1985","Feb"]<- mean(c(Temp["1984","Feb"],Temp["1986","Feb"]))
Temp["1985","Mar"]<- mean(c(Temp["1984","Mar"],Temp["1986","Mar"]))
Temp["2007","Jan"]<- mean(c(Temp["2005","Jan"],Temp["2008","Jan"]))
Temp["2007","Feb"]<- mean(c(Temp["2005","Feb"],Temp["2008","Feb"]))

# Compute annual mean, year 2003 still to be corrected  -------------------

Temp$AnnualMean=rowMeans(Temp[,c(2:13)], na.rm=F)
#Create time series:
AnnualMean.TS=ts(Temp$AnnualMean,start=1969, frequency=1)

# Temperature data from Breitenbach book -------------------------------------

Temp_book= read.csv("Temp_TRAP2_from_Breitenbach_Book.csv",sep=";",row.names=1)
head(Temp_book)
names(Temp_book)
Temp_book$year
#remove first year because incomplete:
Temp_book=Temp_book[-1,]
#Create time series:
Temp_book$year
Temp_book.TS=ts(Temp_book$mean,start=1986, frequency=1)


# Compare with air temperature ----------------------------------------------

AirTemp= read.csv("Breitenbach_AirTemp.csv",sep=",")
head(AirTemp)

AirTemp1969_2010=AirTemp[AirTemp$Year>1968 & AirTemp$Year<2011,]
AirTTS=ts(AirTemp1969_2010$AirT,start=1969, frequency=1)


plot(Temp_book.TS, col="grey", xlim=c(1969,2010), ylim=c(6,10), lwd=6, ylab="Temperature")
lines(AnnualMean.TS, col="red", lwd=2)
lines(AirTTS, col="lightblue", lwd=1)
legend("bottomright", c("Data from Rüdiger", "Data from Breitenbach book","Air T"), col=c("grey", "red", "lightblue"), lty=1, lwd=2)

# change year 2003 and 2006 with values from Breitenbach book ----------------------

Temp["2003","AnnualMean"]<- Temp_book["2003","mean"]
Temp["2006","AnnualMean"]<- Temp_book["2006","mean"]

AnnualMean.TS=ts(Temp$AnnualMean,start=1969, frequency=1)
plot(AirTTS,ylim=c(6,10), col="lightblue")
lines(AnnualMean.TS, lwd=2)
legend("bottomright", c("Water T", "Air T"), col=c("black", "lightblue"), lty=1, lwd=2)

### Note:
### water temperature seems too high in 1984 compared to air temperature.

# Mean annual water temperature - Final table  ----------------------------
WaterTBreitenbach=Temp[, c("year","AnnualMean")]
#export table
write.table(WaterTBreitenbach, "WaterTBreitenbach_final.csv", row.names = FALSE, append = FALSE, col.names = TRUE, sep = ", ", quote = TRUE)
