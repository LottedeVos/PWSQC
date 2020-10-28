##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


# PWS QC - script 1: the construction of the neighbour list of all stations in the PWS dataset

# Input:
# -Text file 'Meta' (located in Folder 'InputFiles') containing the locations in Lon and Lat and unique ID, corresponding to those in the names of the raw PWS observation text files, and column 'radref' with the unique ID of the corresponding radar observation text file, of all available stations with rain observations. 
# -Filtersettings as provided in Inputfiles/FilterSettings.txt.

# Output:
# -A list of length equal to the number of stations in Meta data, of vectors providing the ID's of neighbour stations in the order of ID's in Meta. 

# Notes:
# -The output is written in Folder 'OutputFolder', with the label of that combination of Filtersettings. 
# -When PWSQC_0.R has run the folder "OutputFolder" was created. Make sure this folder exists as the output is written there. 
# -The project() function assumes locations in the Netherlands, and this function needs to be adjusted for other areas in the world where the rijksdriehoekscoordinaten projection is not approriate. 


rm(list=ls())
library(geosphere) # to calculate distance bewteen lat/lon coordinates

workingdirectory <- "..."	# pathway to workingdirectory where the R-scripts are located.
setwd(workingdirectory)
source("InputFiles/Filtersettings.txt")	# obtain 'range', 'nstat', 'nint', 'HIthresA', 'HIthresB', 'compareint', 'rainyint', 'matchint', 'corthres' and 'Filtersettings'


Meta <- read.table("InputFiles/metatableAms.txt", header=T, sep=",")

 # # Construction of neigbourlist for each station # #

neighbourlist <- vector("list", nrow(Meta)) 
for(i in 1:nrow(Meta)){
	dist <- distm(cbind(Meta$lon, Meta$lat), c(Meta$lon[i], Meta$lat[i]), fun = distHaversine)	# make a list of distances to all stations including itself
	neighbourlist[[i]] <- Meta$id[which(dist > 0 & dist <= range)] }	# select the ID's of stations where the distance is smaller than 'range' but larger than zero to avoid being matched with itself

save(neighbourlist, file=paste0("OutputFolder/neighbourlist_Filtersettings",Filtersettings,".RData"))	# save 'neighbourlist' as list in an R Object 
