##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# PWS QC - script 0: the construction of the PWS dataset and Radar reference dataset

# Input:
# -Text file 'Meta' (located in Folder 'InputFiles') containing the locations in Lon and Lat and unique ID, corresponding to those in the names of the raw PWS observation text files, and column 'radref' with the unique ID of the corresponding radar observation text file, of all available stations with rain observations. 
# -Text files with timestamp and rainfall amount columns of rainfall observations by weather stations, where the timestamp indicates the end interval of the observation and the rainfall value indicates the amount of rainfall in mm since last observation. The timeseries may have variable time intervals and can contain large data gaps. In the current example the timecolumn are provided in Unix timestamp. 
# -Text files with timestamp and rainfall amount columns of rainfall observations of the gauge-adjusted radar product, in 5 min resolution at the overlying pixels of the weather stations. In the current example the timecolumn are provided in YYYYmmddHHMMSS in GMT. 

# Output: 
# -A matrix 'Ndataset' with columns of 5 min timeseries between 'starttime' and 'endtime' containing the most recent measurements from the Netatmo dataset at each timestamp and 'NA' if no observations was available. The order of the columns corresponds with the order of the IDs in the 'Meta' table. This is saved as R Object. 
# -A matrix 'Rdataset' with columns of 5 min timeseries between 'starttime' and 'endtime' containing the measurements from the radar reference dataset since last 5 min timestamp and 'NA' if no observations was available. The order of the columns corresponds with the order of the IDs in the 'Meta' table. This is saved as R Object. 

# Notes: 
# -The text files that were used as input have already been produced. The meta file includes meta-data of all stations. It doesn't matter if a station contains no observations for the study period. The radar timeseries have been produced based on the Lon and Lat of the weather stations in the meta file. 
# -Adjust the pathway according to the location of the files on your computer. 
# -Because 'Ndataset' and 'Rdataset' are saved as R-objects in directory 'OutputFolder', this script only needs to be run once. 
# -Note that the SO-filter and biascorrectionfactor estimation require a start-up period, which depending on the Filtersettings, rain occurrence and data availability may take up to a month. For this reason we run the filters for a 13 month period and exclude the first month in the analysis of the results. 
# -The Meta file we use in our example only includes stations with rain observations.
# -The raw Radar observations were stored in combined txt-files containing observations of 10 different station ID's in the example case. For that reason, obtaining the relevant observations from the text files is slightly more complex than simply calling the file by it's ID. This should be adjusted accordingly when applying the scripts on a different dataset. Also adjust as needed for name formats and format of input (e.g. "xx.rain.historic.csv" and "radar_5min_Netatmo_xx", seperated by ";" or ",", names of colums, etc.). 



rm(list=ls())
require(data.table)	# needed for setDT()

workingdirectory <- "..."	# pathway to workingdirectory where the R-scripts are located.
NfilesFolder <- "..." # pathway to the folder where the weather station timeseries files are located. 
RfilesFolder <- "..." # pathway to the folder where the radar reference timeseries files are located. 
setwd(workingdirectory)
dir.create("OutputFolder")


 # # period and station information # #

starttime <- strptime(20160501, format="%Y%m%d", tz="GMT")	# change if needed
endtime <- strptime(20170601, format="%Y%m%d", tz="GMT")	# change if needed
Time <- seq(starttime, endtime, by="5 min")[-1]

Meta <- read.table("InputFiles/metatableAms.txt", header=T, sep=",")


 # # construct Ndataset # #

for(i in 1:nrow(Meta)){
	print(paste("Ndataset construction progress:", i, "out of", nrow(Meta)))
	N <- read.table(paste0(NfilesFolder,Meta$id[i],".rain.historic.csv"), header=T, sep=";")
	if(nrow(N) < 2){ Ncomplete <- rep(NA, times=length(Time))	# if no observations are available, make a sequence of NA-values
		}else{
	N$timestamp <- as.POSIXct(N$timestamp, origin="1970-01-01", tz="GMT")	
	if(range(N$timestamp)[1] > endtime | range(N$timestamp)[2] < starttime){Ncomplete <- rep(NA, times=length(Time))	# if 'N' does not have data within timerange, make a sequence of NA-values
		}else{
	if(length(which(is.na(N$rain_mm)))>0){N <- N[-which(is.na(N$rain_mm)),]}	# exclude observations in 'N' that are NA
	N <- unique(N[,1:2])	# exclude duplicated rows
	N <- N[order(N$timestamp),] 	# make sure 'N' is chronological 
	N <- N[which(N$timestamp >= starttime & N$timestamp < endtime),]  # select only observations within the study period	

	Ntimeround <- strptime(cut(c(starttime, N$timestamp), "5 min")[-1], format="%Y-%m-%d %H:%M", tz="GMT")	+ 5*60	# because we apply the filter every 5 min, 'starttime' is included to force start at round 5 min value; 5 minutes are added because cut() rounds down
	Nagg <- setDT(as.data.frame(N$rain_mm))[,lapply(.SD,sum), by=.(Ntimeround)]	# take the sum of observations within the same 5 min interval
	names(Nagg) <- c("Time", "Rain")    # name-giving is required to merge with the 'Time' column in next line
	Ncomplete <- merge(as.data.frame(Time), Nagg, all=T)$Rain	# make sure that there is an observation for each timestep
	if(length(Ncomplete) != length(Time)){print(paste0("Ncomplete length ERROR:",i))}   # give a warning if for some reason 'Ncomplete' is not equally long to 'Time'
	}} #end else loops
   if(exists("Ndataset")==F){ Ndataset <- Ncomplete
	}else{ Ndataset <- cbind(Ndataset, Ncomplete) } # add 'Ncomplete' to 'Ndataset' as new column
} # end i-loop 	

save(Ndataset, file=paste0("OutputFolder/Ndataset.RData"))	# save 'Ndataset' as matrix in an R Object 


 # # construct Rdataset # #

for(i in 1:nrow(Meta)){
	print(paste("Rdataset construction progress:", i, "out of", nrow(Meta)))
   	R <- read.table(paste0(RfilesFolder,"radar_5min_Netatmo_",ifelse(nchar(ceiling(Meta$radref[i]/10)) == 1, paste0(0, ceiling(Meta$radref[i]/10)), ceiling(Meta$radref[i]/10)),".csv"), header=F, sep=",")
   	R <- R[which(R[,1] == Meta$radref[i]),2:3]
   	names(R) <- c("Time", "Rain")
   	R$Time <- as.POSIXct(format(R$Time, scientific=F), format="%Y%m%d%H%M%S", tz="GMT")
   	R <- R[which((R$Time > starttime)&(R$Time <= endtime)),]  # select only observations within the study period		

   	if(exists("Rdataset")==F){ Rdataset <- R$Rain
		}else{ Rdataset <- cbind(Rdataset, R$Rain) }
} # end i-loop 	

save(Rdataset, file=paste0("OutputFolder/Rdataset.RData"))	# save 'Rdataset' as matrix in an R Object 

