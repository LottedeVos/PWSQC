##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# PWS QC - script 3: the attribution of SO flags and bias correction factors to the Ndataset

# Input:
# -'Ndataset' as constructed in PWSQC_0.R.
# -'neighbourlist' as constructed in PWSQC_1.R.
# -'FZ_flags' as constructed in PWSQC_2.R.
# -'HI_flags' as constructed in PWSQC_2.R.
# -Text file 'Meta' (located in Folder 'InputFiles') containing the locations in Lon and Lat and unique ID, corresponding to those in the names of the raw PWS observation text files, and column 'radref' with the unique ID of the corresponding radar observation text file, of all available stations with rain observations. 
# -Text file 'Filtersettings' (located in Folder 'InputFiles')

# Output: 
# -A matrix 'SO_flags' with equal size as 'Ndataset' that consists of '0', '1' or '-1' to indicate the attribution of Station Outlier flags. The order of the columns corresponds with the order of the ID's in the meta table. This is saved as R Object. 
# -A matrix 'biascorrectiontable' with equal size as 'Ndataset' that consists of any value needed to compensate the bias for that station at that time. The order of the columns corresponds with the order of the ID's in the meta table. This is saved as R Object. 

# Notes: 
# -Warning: this analysis may take over two days runtime for the example year dataset, as the past period is evaluated at each timestep.  
# -The construction of 'Ndataset' and 'Rdataset' does not rely on any of the filter settings. They do depend on the period. Make sure the chosen 'starttime' and 'endtime' correspond with those in PWSQC_0.R.
# -All stations start out with the 'defaultbiascorrection' provided in the Filtersettings.txt file, and adjust accordingly after a lead-time of ~1 month, depending settings, rain occurrence and data availability. 
# -The Station Outlier is determined based on the intervals that were not flagged as HI and/or FZ.
#  0 = No flag
#  1 = A flag is attributed
# -1 = Not enough information to determine flag 


rm(list=ls())

workingdirectory <- "..."	# pathway to workingdirectory where the R-scripts are located
setwd(workingdirectory)
source("InputFiles/Filtersettings.txt")	# obtain 'range', 'nstat', 'nint', 'HIthresA', 'HIthresB', 'compareint', 'rainyint', 'matchint', 'corthres' and 'Filtersettings'
load("OutputFolder/Ndataset.RData")
load(paste0("OutputFolder/FZ_flags_Filtersettings",Filtersettings,".RData"))
load(paste0("OutputFolder/HI_flags_Filtersettings",Filtersettings,".RData"))
load(paste0("OutputFolder/neighbourlist_Filtersettings",Filtersettings,".RData"))


 # # station information # #

Meta <- read.table("InputFiles/metatableAms.txt", header=T, sep=",")

Ndataset2 <- Ndataset * defaultbiascorrection	# the multiplicationfactor does not impact SO in any way 
Ndataset2[which((HI_flags == 1)|(FZ_flags == 1))] <- NA


 # # Construct SO flags and biascorrectiontable # #

for(i in 1:nrow(Meta)){
	Nint <- Ndataset2[,i]

   if((length(neighbourlist[[i]]) < nstat)|(length(which(is.na(Nint)==F)) < 1)){SOflag <- rep(-1, times=length(Nint))	# if there are not enough stations nearby or no observations in 'Nint', all intervals get flagged as -1
	}else{

	Nintrain <- rep(0, length=length(Nint))
	Nintrain[which(Nint > 0)] <- 1
	Nintraincum <- cumsum(Nintrain)	# cumulative intervals with nonzero rainmeasurements
	comparestartrowA <- match((Nintraincum-rainyint+1), Nintraincum)-1	# row from which the window should start to have at least 'rainyint' rainy intervals. match the first value that belongs to a number of rainy intervals equal to 'rainyint'-1 (double minus becomes +), and detract 1 after. This makes sure that if multiple rows have same amount, the last one (nearest to matching interval) is chosen. Note that this may result in a rownumber of 0 (1 minus 1), which should be replaced by NA (see next line), as the 0th row does not exist. 
	comparestartrowA[which(comparestartrowA == 0)] <- NA
	comparestartrowB <- c(rep(NA, times=compareint-1), 1:(length(Nint)-compareint+1))	# row from which the window should start to have 'compareint' number of intervals
	comparestartrow <- ifelse(is.na(comparestartrowB), NA, ifelse((comparestartrowA < comparestartrowB), comparestartrowA, comparestartrowB))	# choose either 'compareint' steps before, or where 'rainyint' was reached 	

   NeighbourVal <- Ndataset2[,which(Meta$id %in% neighbourlist[[i]])]	
   NeighbourVal[which(is.na(Nint)),] <- NA		# replace rows where 'Nint' is NA with NA, this is needed for checking 'matchint' later

   cortable <- biastable <-  matrix(NA, ncol=ncol(NeighbourVal), nrow=nrow(NeighbourVal)) # table of size 'NeighbourVal' to fill in the correlation values

   	for(t in 1:length(Nint)){
	print(paste("SO filter: station", i, "out of", nrow(Meta), " -  t",t, "out of", length(Nint))) 
		if(is.na(comparestartrow[t])){next}
		NeighbourValselec <- NeighbourVal[comparestartrow[t]:t,]
		columnselec <- which(apply(NeighbourValselec, 2, function(x) length(which(is.na(x)==F))) > matchint)	# find columns that have at least 'matchint' overlapping intervals
		if(length(columnselec) < nstat){next}
		cortable[t,columnselec] <- apply(NeighbourValselec[,columnselec], 2, function(x) cor(x, Nint[comparestartrow[t]:t], use='complete.obs')) # determine the correlations with all neighbour stations over the past period that was chosen, this can yield a NA when the comparing stations measures zeroes only
		biastable[t,columnselec] <- apply(NeighbourValselec[,columnselec], 2, function(x) mean(Nint[comparestartrow[t]:t]/defaultbiascorrection - x, na.rm=T)/mean(x, na.rm=T) ) # calculate the relative bias in the mean of the raw observations with all biascorrected neighbour stations over the past period between 'comparestartrow[t]' and 't'. To do so, 'Nint', which is based on 'Ndataset' times 'defaultbiascorrection', should be divided by 'defaultbiascorrection'.
   	} # end of r-loop  

   SOflag <- rep(0, times=length(Nint))
   SOflag[which(apply(cortable, 1, function(x) median(x, na.rm=T)) < corthres)] <- 1	# rows where the median of all neighbour correlations was below 'corthres' are flagged as Station Outlier
   SOflag[which(apply(cortable, 1, function(x) length(which(is.na(x)==F))) < nstat)] <- -1	# 'cortable' will have more or equal amounts of NA values than 'biastable'. Where 'SOflag' is 0, enough stations were included in the calculation to consider that row in 'biastable'.
	} # end of ifelse loop 
   if(exists("SO_flags")==F){ SO_flags <- SOflag
	}else{ SO_flags <- cbind(SO_flags, SOflag) }

   biascorrectiontimeline <- rep(defaultbiascorrection, times=length(Nint)) # start out with the 'defaultbiascorrection' for all stations at all time intervals
   if(length(which(SOflag == 0)) > 0){	
   biasmed <- apply(biastable, 1, function(x) median(x, na.rm=T))	# the median of the bias between raw timeserie 'Nint' with neighbouring stations that are multiplied with the 'defaultbiascorrection'

   for(brow in which(SOflag == 0)){
	biasprev <- biascorrectiontimeline[brow]
	biasnew <- 1 / (1+biasmed[brow])
	if( abs(log(biasnew / biasprev)) > log(1+biasthres) ){	# means that if [1/(1+biasthres) > BCFnew/BCFprev > 1+biasthres], change it for the remainder of the timeline   
		biascorrectiontimeline[(brow+1):length(biascorrectiontimeline)] <- biasnew }
   } # end of brow-loop
   } # end of if-loop	

   if(exists("biascorrectiontable")==F){ biascorrectiontable <- biascorrectiontimeline
	}else{ biascorrectiontable <- cbind(biascorrectiontable, biascorrectiontimeline) }

} # end of i loop

save(SO_flags, file=paste0("OutputFolder/SO_flags_Filtersettings",Filtersettings,".RData"))	# save 'SO_flags' as matrix in an R Object 
save(biascorrectiontable, file=paste0("OutputFolder/biascorrectiontable_Filtersettings",Filtersettings,".RData"))	# save 'biascorrectiontable' as matrix in an R Object 

