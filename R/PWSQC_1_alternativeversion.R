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
# With respect to this last point:
# Initially the code made a projection of the lon and lat locations (in rijksdriehoekscoordinaten) in order to calculate distances between stations in meter. 
# The projection that was used is appropriate for the Netherlands, but not for other areas in the world. 
# These lines in the code had been replaced with a function from the geosphere package, which calculates distance between stations directly from their lon and lat locations. 
# The function (distHaversine) assumes a spherical earth, ignoring ellipsoidal effects. The new code is applicable on PWS networks anywhere in the world.
# Courtesy of Dr. Maarten Reyniers and Eva Beele.
# Much faster solution is to use the sf package. This seems always to give -0.1117411% shorter distances, but seems 15 times faster, and has therefore been implemented below.
# In order to save time, the number of neighbour stations is limited by a maximimum "MaxNrStations", if more nearby stations are available only the closest neighbours are considered. 
# Courtesy of Dr. Aart Overeem and Niek van Andel.

rm(list=ls())
library(geosphere,lib.loc=path) # to calculate distance between lat/lon coordinates
library(sf,lib.loc=path)
library(s2,lib.loc=path)

workingdirectory <- "..."	# pathway to workingdirectory where the R-scripts are located.
setwd(workingdirectory)
source("InputFiles/Filtersettings.txt")	# obtain 'range', 'nstat', 'nint', 'HIthresA', 'HIthresB', 'compareint', 'rainyint', 'matchint', 'corthres' and 'Filtersettings'

MaxNrStations <- "..." # the maximum number of neighbours within the range to be considered, e.g. 20

Meta <- read.table("InputFiles/metatableAms.txt", header=T, sep=",")


 # # Construction of neigbourlist for each station # #

neighbourlist <- vector("list", nrow(MetaTotal)) 

pnts <- data.frame(id = Meta$id, lon = Meta$lon, lat = Meta$lat)
pnts_sf <- st_as_sf(pnts, crs = 4326L, coords = c("lon", "lat")) # crs = 4326L is EPSG:4326. WGS 84 -- WGS84 - World Geodetic System 1984, used in GPS  


LonLat <- cbind(Meta$lon, Meta$lat)
for(i in 1:nrow(Meta)){
        dist <- as.numeric(st_distance(pnts_sf, pnts_sf[i,])) # make a list of distances to all stations including itself
        distSelected <- dist[which(dist > 0 & dist <= range)]
        # Determine maximum range of MaxNrStations nearest stations in order to select the MaxNrStations nearest stations:
        if (length(distSelected) > MaxNrStations)
        {
        	range_max <- max(sort(distSelected)[1:MaxNrStations])
        }
        else
        {
	        range_max <- range
        }
	neighbourlist[[i]] <- Meta$id[which(dist > 0 & dist <= range_max)] 	# select the ID's of stations where the distance is smaller than 'range' but larger than zero to avoid being matched with itself
}

save(neighbourlist, file=paste0("OutputFolder/neighbourlist_Filtersettings",Filtersettings,".RData"))	# save 'neighbourlist' as list in an R Object 




