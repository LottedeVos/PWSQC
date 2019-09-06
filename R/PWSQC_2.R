##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# PWS QC - script 2: the attribution of FZ and HI flags to the Ndataset

# Input:
# -'Ndataset' as constructed in PWSQC_0.R.
# -'neighbourlist' as constructed in PWSQC_1.R.
# -Text file 'Meta' (located in Folder 'InputFiles') containing the locations in Lon and Lat and unique ID, corresponding to those in the names of the raw PWS observation text files, and column 'radref' with the unique ID of the corresponding radar observation text file, of all available stations with rain observations. 
# -Text file 'Filtersettings' (located in Folder 'InputFiles')

# Output: 
# -A matrix 'FZ_flags' with equal size as 'Ndataset' that consists of '0', '1' or '-1' to indicate the attribution of Faulty Zero flags. The order of the columns corresponds with the order of the IDs in the Meta table. This is saved as R Object. 
# -A matrix 'HI_flags' with equal size as 'Ndataset' that consists of '0', '1' or '-1' to indicate the attribution of High Influx flags. The order of the columns corresponds with the order of the IDs in the Meta table. This is saved as R Object. 

# Notes: 
# -The construction of 'Ndataset' and 'Rdataset' does not rely on any of the filter settings. They do depend on the period. Make sure the chosen 'starttime' and 'endtime' correspond with those in PWSQC_0.R.
# -Allthough the HI-filter and FZ-filter are seperate modules, these flags are attributed at the same time for time efficiency.
#  0 = No flag
#  1 = A flag is attributed
# -1 = Not enough information to determine flag 


rm(list=ls())

workingdirectory <- "..."	# pathway to workingdirectory where the R-scripts are located
setwd(workingdirectory)
source("InputFiles/Filtersettings.txt")	# obtain 'range', 'nstat', 'nint', 'HIthresA', 'HIthresB', 'compareint', 'rainyint', 'matchint', 'corthres' and 'Filtersettings'
load("OutputFolder/Ndataset.RData")
load(paste0("OutputFolder/neighbourlist_Filtersettings",Filtersettings,".RData"))


 # # station information # #

Meta <- read.table("InputFiles/metatableAms.txt", header=T, sep=",")


 # # Construct HI and FZ flags # #

for(i in 1:nrow(Meta)){
   print(paste("FZ_flags and HI_flags construction progress:", i, "out of", nrow(Meta)))

	Nint <- Ndataset[,i]
	if((length(which(is.na(Nint)==F)) < 1) | (length(neighbourlist[[i]]) < nstat)){	# if the 'Nint' column consist of no observations or there are too few neighbours, make a sequence of -1 values
		  HIflag <- FZflag <- rep(-1, times=length(Nint))
   		  if(exists("HI_flags")==F){ HI_flags <- HIflag
			}else{ HI_flags <- cbind(HI_flags, HIflag) }
 		  if(exists("FZ_flags")==F){ FZ_flags <- FZflag
			}else{ FZ_flags <- cbind(FZ_flags, FZflag) }
	}else{

   NeighbourVal <- Ndataset[,which(Meta$id %in% neighbourlist[[i]])]	# take a subset of 'Ndatasset' with only the columns corresponding with the ID's of the neighbouring stations

   Ref <- rep(NA, times=length(Nint))
   Number_of_measurements <- apply(NeighbourVal, 1, function(x) length(which(is.na(x)==F)))	# count the number of neighbours with measurements at each interval
   Med <- apply(NeighbourVal, 1, median, na.rm=T)	# take the median of all neighbour values

   # # # HI-filter:
   HIflag <- rep(0, times=length(Nint))
   HIflag[which(((Nint > HIthresB) & (Med < HIthresA)) | ((Med >= HIthresA) & (Nint > (HIthresB*Med/HIthresA))))] <- 1  # if thresholds are exceeded, the HI flag becomes 1
   HIflag[which(Number_of_measurements < nstat)] <- -1	# if less than 'nstat' neighbours supply observations, the HI flag becomes equal to -1

   if(exists("HI_flags")==F){ HI_flags <- HIflag
	}else{ HI_flags <- cbind(HI_flags, HIflag) }

   # # # FZ-filter:
   Ref[which(Med == 0)] <- 0
   Ref[which(Med >  0)] <- 1
   Ref[which(Number_of_measurements < nstat)] <- NA	# make binary reference with 1 for rainy periods and 0 for dry periods (based on median of neighbour observations and NA when not enough info is available

   Nwd <- Nint
   Nwd[which(Nint > 0)] <- 1		# make binary timeseries of station observations where 1 stands for wet and 0 for dry observations
   runs <- rle(Nwd)
   rownr <- cumsum(runs$lengths)	
   endrow <- rownr[ which(runs$lengths > nint & runs$values==0) ]
   startrow <- endrow - runs$lengths[ which(runs$lengths > nint  & runs$values==0) ] + 1	# 'endrow' and 'startrow' indicate the boundaries of dry periods as measured by the station

   FZflag <- rep(0, times=length(Nint))	
   if(length(endrow) > 0){
   for(r in 1:length(endrow)){
   	if(length( which( (Ref[startrow[r] : endrow[r]]) == 1) ) > nint ){	# in case at least 'nint' intervals in 'Ref' are wet where 'Nint' is dry

	runs2 <- rle(Ref[startrow[r] : endrow[r]])	# check if the 'nint' wet intervals in 'Ref' are consecutive. 
   	rownr2 <- cumsum(runs2$lengths)	
   	endrow2 <- rownr2[ which(runs2$lengths > nint & runs2$values==1) ]
   	startrow2 <- endrow2 - runs2$lengths[ which(runs2$lengths > nint  & runs2$values==1) ] + 1	

	if(length(startrow2) > 0){
	FZstartrow <- startrow[r] + startrow2[1] - 1 + nint	# the interval in 'Nint' where the previous 'nint' intervals were dry in 'Nint' and wet in the median

   	FZflag[FZstartrow : endrow[r]] <- 1	# from this interval up to the end of the dry period is flagged as Faulty Zero
	
	m <- 1
	while((is.na(Nwd[endrow[r] + m])|(Nwd[endrow[r] + m] == 0)) & ((endrow[r]+m) <= length(Nwd)) ){ # if subsequent values in 'Nwd' are NA or 0, continue flagging until 'Nwd' becomes 1 or the end of 'Nwd' is reached
	 FZflag[endrow[r]+m] <- 1	# once a period is labeled as Faulty Zero, flagging continues until there is a rain measurement (NA values are ignored)
	 m <- m+1
	} # end while-loop

	}} # end if loops
   } # end r-loop
   } # end if loop (endrow)

   FZflag[which(Number_of_measurements < nstat)] <- -1	# if too few neighbours have observations the FZ can not be attributed
   if(exists("FZ_flags")==F){ FZ_flags <- FZflag
	}else{ FZ_flags <- cbind(FZ_flags, FZflag) }

	} #end of ifelse loop

} # end of i-loop

save(FZ_flags, file=paste0("OutputFolder/FZ_flags_Filtersettings",Filtersettings,".RData"))	# save 'FZ_flags' as matrix in an R Object 
save(HI_flags, file=paste0("OutputFolder/HI_flags_Filtersettings",Filtersettings,".RData"))	# save 'HI_flags' as matrix in an R Object 





