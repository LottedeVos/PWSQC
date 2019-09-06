##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# PWS QC - script 4: validating the results of the QC using gauge-adjusted radar reference

# Input:
# -'Ndataset' as constructed in PWSQC_0.R.
# -'RNdataset' as constructed in PWSQC_0.R.
# -'FZ_flags' as constructed in PWSQC_2.R.
# -'HI_flags' as constructed in PWSQC_2.R.
# -'SO_flags' as constructed in PWSQC_3.R.
# -'biascorrectiontable' as constructed in PWSQC_3.R.
# -Text file 'Meta' (located in Folder 'InputFiles') containing the locations in Lon and Lat and unique ID, corresponding to those in the names of the raw PWS observation text files, and column 'radref' with the unique ID of the corresponding radar observation text file, of all available stations with rain observations. 
# -Text file 'Filtersettings' (located in Folder 'InputFiles')

# Output: 
# -A doublemass plot of PWS cumulative rainfall vs radar reference cumulative rainfall before and after QC is applied, indicating the flagged intervals by color. The fraction of remaining intervals and the correlation with radar reference values are calculated for filtered timelines and included in a inset boxplots figure. 

# Notes: 
# -In the example a period of 13 months is evaluated but only validate for the last 12 month period, as we consider the first month as lead-up time. 
# -Here we apply the bias correction, however, one may choose to adjust the scripts in PWSQC_3.R to not generate a biascorrectionfactor but a flag-warning when the bias changes in time. 
# -One can make the results by only excluding the flagged intervals (flag = 1) or excluding all intervals that are not certain to be okay (flag = 1 and -1) by making 'FiltersettingStrict' either 'no' or 'yes'. 
#  0 = No flag
#  1 = A flag is attributed
# -1 = Not enough information to determine flag 
# in red: the FZ flags
# in orange: the HI flags
# in green: the SO flags
# overall adjustment with biascorrection factor
# the CV and cor values are calculated for filtered timelines and included in a inset boxplots figure. 


rm(list=ls())
require(TeachingDemos)	# needed for subplot()

workingdirectory <- "..."	# pathway to workingdirectory where the R-scripts are located
setwd(workingdirectory)
source("InputFiles/Filtersettings.txt")	# obtain 'range', 'nstat', 'nint', 'HIthresA', 'HIthresB', 'compareint', 'rainyint', 'matchint', 'corthres' and 'Filtersettings'
load("OutputFolder/Ndataset.RData")
load("OutputFolder/Rdataset.RData")
load(paste0("OutputFolder/FZ_flags_Filtersettings",Filtersettings,".RData"))
load(paste0("OutputFolder/HI_flags_Filtersettings",Filtersettings,".RData"))
load(paste0("OutputFolder/SO_flags_Filtersettings",Filtersettings,".RData"))
load(paste0("OutputFolder/biascorrectiontable_Filtersettings",Filtersettings,".RData"))


FiltersettingStrict <- "yes"		# 'yes' means that only flag == 0 intervals are included. 'no' means that flag == 0 and flag == -1 are included 
lim <- c(0,1000)  # if needed change upper boundary to maximum cumulative rainfall over the selected period

 # # period and station information # #
starttime_with_leadtime <- strptime(20160501, format="%Y%m%d", tz="GMT")
starttime <- strptime(20160601, format="%Y%m%d", tz="GMT")
endtime <- strptime(20170601, format="%Y%m%d", tz="GMT")

Meta <- read.table("InputFiles/metatableAms.txt", header=T, sep=",")

	# leave out the first month lead-time in the figures:
	Time <- seq(starttime_with_leadtime, endtime, by="5 min")[-1]
	Ndataset <- Ndataset[which(Time > starttime)[1] : which(Time == endtime),]
	Rdataset <- Rdataset[which(Time > starttime)[1] : which(Time == endtime),]
	FZ_flags <- FZ_flags[which(Time > starttime)[1] : which(Time == endtime),]
	HI_flags <- HI_flags[which(Time > starttime)[1] : which(Time == endtime),]
	SO_flags <- SO_flags[which(Time > starttime)[1] : which(Time == endtime),]
	biascorrectiontable <- biascorrectiontable[which(Time > starttime)[1] : which(Time == endtime),]

png(paste0("OutputFolder/doublemassplot_Ams_HI_FZ_SO_bias_filterstrict",FiltersettingStrict,"_filtersettings",Filtersettings,".png"), width=1080, height=540)

par(mar=c(5.1,5.1,4.1,1.1))
layout(matrix(c(1,2),1,2,TRUE))	# needed to make subplot in right window
plot(NA, NA, xlab="Gauge-adjusted radar cumulative rain (mm)", ylab="PWS cumulative rain (mm)", xlim=lim, ylim=lim, main="Raw", cex.lab=2, cex.axis=1.5, cex.main=1.8)
abline(0, 1, col="grey", lty=2, lwd=2)

for(i in 1:nrow(Meta)){
	R <- Rdataset[,i]
	N <- Ndataset[,i]
	selec <- which(is.na(R) | is.na(N))
		FZflag <- FZ_flags[-selec,i]
		HIflag <- HI_flags[-selec,i]
		SOflag <- SO_flags[-selec,i]

	N <- N[-selec]
	R <- R[-selec]

	Raccum <- cumsum(R)
	Naccum <- cumsum(N)

	print(paste("i =", i, "; ID =", Meta$id[i]))

	lines(Raccum, Naccum, lwd=0.5, col="blue")
   	points(Raccum[which(FZflag == 1)], Naccum[which(FZflag == 1)], col="red", pch=1, cex=0.5)
	segments(x0=Raccum[which(HIflag == 1)], x1=Raccum[which(HIflag == 1)-1], y0=Naccum[which(HIflag == 1)], y1=Naccum[which(HIflag == 1)-1], col="orange", lwd=3)
	points(Raccum[which(SOflag == 1)], Naccum[which(SOflag == 1)], col="green", pch=1, cex=0.5)	

}#end i-loop
	legend(max(lim)*0.1, max(lim), c("FZ-flag", "HI-flag", "SO-flag"), pch=16, col=c("red", "orange", "green"), cex=1.2)

corlist <- CVlist <- c()
fracremain <- rep(0, times=nrow(Meta))

plot(NA, NA, xlab="Gauge-adjusted radar cumulative rain (mm)", ylab="PWS cumulative rain (mm)", xlim=lim, ylim=lim, main="Filtered", cex.lab=2, cex.axis=1.5, cex.main=1.8)
abline(0, 1, col="grey", lty=2, lwd=2)

for(i in 1:nrow(Meta)){
	R <- Rdataset[,i]
	N <- Ndataset[,i]
	FZflag <- FZ_flags[,i]
	HIflag <- HI_flags[,i]
	SOflag <- SO_flags[,i]
	biascorrection <- biascorrectiontable[,i]

	if(FiltersettingStrict == 'yes'){filtered <- which(is.na(R) == F & is.na(N) == F & FZflag == 0 & HIflag == 0 & SOflag == 0)}
	if(FiltersettingStrict == 'no' ){filtered <- which(is.na(R) == F & is.na(N) == F & FZflag != -1 & HIflag != -1 & SOflag != -1)}

	N <- (N*biascorrection)[filtered]
	R <- R[filtered]

		if(length(N) > 1){
		corlist[i] <- cor(N, R)
		CVlist[i] <- sd(N - R) / mean(R)
		fracremain[i] <- length(filtered) / (nrow(Ndataset) - length(which(is.na(R) | is.na(N))))	# fraction filtered over fraction that was not NA
		}

	Raccum <- cumsum(R)
	Naccum <- cumsum(N)

	print(paste("i =", i, "; ID =", Meta$id[i]))

	lines(Raccum, Naccum, lwd=0.5, col="blue")
}#end i-loop

rect(max(lim)*0.05, max(lim)*0.6,  max(lim)*0.4, max(lim), col = "white", border=NA)	# make white background for subplot
subplot( boxplot(corlist, fracremain, names=c("correlation","fraction"), cex.axis=1, col=c("lightslateblue", "lightblue"), ylim=c(0,1), yaxs="i"), x=c(max(lim)*0.05,max(lim)*0.4), y=c(max(lim)*0.6,max(lim)))

dev.off()


# # calculate validation metrics # #
	nodata <- which(is.na(as.vector(Rdataset)) | is.na(as.vector(Ndataset))) # as intervals without measurements can have flag = 0 or -1 depending on whether they will have measurements later in the timeline
	frac1 <- 100
	frac2 <- round(100* length(which(FZ_flags[-nodata] != 1 & HI_flags[-nodata] != 1 & SO_flags[-nodata] != 1)) / length(Ndataset[-nodata]), digits=3)	# fraction on non-NA data that was not flagged (Flex)	
	frac3 <- round(100* length(which(FZ_flags[-nodata] == 0 & HI_flags[-nodata] == 0 & SO_flags[-nodata] == 0)) / length(Ndataset[-nodata]), digits=3)	# fraction on non-NA data that was not flagged (Strict)	
	
  	bias1 <- round( mean(as.vector(Ndataset)[-nodata] - as.vector(Rdataset)[-nodata]) / mean(as.vector(Rdataset)[-nodata]), digits=3)	# relative bias in the mean
  	CV1   <- round( sd(as.vector(Ndataset)[-nodata] - as.vector(Rdataset)[-nodata]) / mean(as.vector(Rdataset)[-nodata]), digits=3)		# coefficient of variation
	cor1  <- round( cor(as.vector(Ndataset)[-nodata], as.vector(Rdataset)[-nodata], use='complete.obs'), digits=3)						# pearson correlation

  Ndataset2 <- Ndataset*biascorrectiontable
  Ndataset2[which(FZ_flags == 1 | HI_flags == 1 | SO_flags == 1)] <- NA	# exclude all intervals where a flag is attributed
	nodata2 <- which(is.na(as.vector(Rdataset)) | is.na(as.vector(Ndataset2))) 
  	bias2 <- round( mean(as.vector(Ndataset2)[-nodata2] - as.vector(Rdataset)[-nodata2]) / mean(as.vector(Rdataset)[-nodata2]), digits=3)	# relative bias in the mean
  	CV2   <- round( sd(as.vector(Ndataset2)[-nodata2] - as.vector(Rdataset)[-nodata2]) / mean(as.vector(Rdataset)[-nodata2]), digits=3)		# coefficient of variation
  	cor2  <- round( cor(as.vector(Ndataset2)[-nodata2], as.vector(Rdataset)[-nodata2], use='complete.obs'), digits=3)						# pearson correlation

	Ndataset3 <- Ndataset*biascorrectiontable
	Ndataset3[which(FZ_flags != 0 | HI_flags != 0 | SO_flags != 0)] <- NA	# exclude all intervals where the flag is not zero
	nodata3 <- which(is.na(as.vector(Rdataset)) | is.na(as.vector(Ndataset3))) 
  	bias3 <- round( mean(as.vector(Ndataset3)[-nodata3] - as.vector(Rdataset)[-nodata3]) / mean(as.vector(Rdataset)[-nodata3]), digits=3)	# relative bias in the mean
  	CV3   <- round( sd(as.vector(Ndataset3)[-nodata3] - as.vector(Rdataset)[-nodata3]) / mean(as.vector(Rdataset)[-nodata3]), digits=3)		# coefficient of variation
  	cor3  <- round( cor(as.vector(Ndataset3)[-nodata3], as.vector(Rdataset)[-nodata3], use='complete.obs'), digits=3)						# pearson correlation

	
	print(paste("Unfiltered comparison between Ndataset and Rdataset:    bias =", bias1, "   CV =", CV1, "   cor =", cor1, "   frac =", frac1, "%"))
	print(paste("Filtered Flex comparison between Ndataset and Rdataset:    bias =", bias2, "   CV =", CV2, "   cor =", cor2, "   frac =", frac2, "%"))
	print(paste("Filtered Strict comparison between Ndataset and Rdataset:    bias =", bias3, "   CV =", CV3, "   cor =", cor3, "   frac =", frac3, "%"))





