##| creatieve commons attribution share alike | cc-by-sa 4.0 | Lotte de Vos |##

##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

## 	 These scripts are related to the manuscript 'Quality Control for Crowdsourced Personal Weather Stations to Enable Operational Rainfall Monitoring' by De Vos et al. (2019), see https://doi.org/10.1029/2019GL083731


 # # # # # # # # # # #
 # About the scripts #
 # # # # # # # # # # #

 
# These scripts:
# 1. make a dataset 'Ndataset' with columns of 5 min timeseries between 'starttime' and 'endtime' containing the most recent measurements from the PWS dataset at each timestamp, NA if no observations were available
# 2. determine the lists of neighbouring stations for each station in the dataset based on the predetermined range
# 3. construct two datasets of size 'Ndataset' with all flags* of respectively Faulty Zeroes (FZ) and High Influx (HI), constructed on raw data  
# 4. construct dataset of size 'Ndataset' with Station Outlier (SO) flags* determined on dataset without the intervals where HI and/or FZ flag were 1  
# 5. apply biascorrection and adjusts it when the bias of a stations changes over time
# 6. construct the 5 min timeseries of the gauge adjusted radar data at each station location used as reference
# * [0 is no flag, 1 is flag and -1 is no flag could be given]

# Output is:
# a. Ndataset -> columns of raw Netatmo timeseries where number of columns is the number of stations and number of rows is the number of 5 min intervals between 'tstart' and 'tend' where each timestamp indicates the end of the 5 min interval. 
# b. FZ_flags -> matrix of size Ndataset with FZ flags (-1, 0 or 1)
# c. HI_flags -> matrix of size Ndataset with HI flags (-1, 0 or 1)
# d. SO_flags -> matrix of size Ndataset with SO flags (-1, 0 or 1)
# e. biascorrectiontable -> matrix of size Ndataset with all biascorrections for all timesteps
# f. Rdataset -> columns of 5 min gauge adjusted radar timeseries corresponding with Ndataset



 # # # # # # # # # # # # # #
 # How to use the scripts  #
 # # # # # # # # # # # # # #
 
# More details on the code are provided within the R-scripts. While running the scripts, the progress is tracked by printed statements. Note that especially step 3 is very timeconsuming. All (intermediate) results are saved in OutputFolder. 


  # # preparation # #

## step 0: rewrite timeseries of variable sizes and possible large data gaps into consistent 5 min timeseries
# Run PWSQC_0.R
# notes: File needs to be adjusted according to the pathway to the files. In the current form it assumes files of radar timeseries without data gaps and at 5 min timeseries. Running this script does not have to be repeated when the parameter settings are changed. The name of the Meta file in this example is assumed to be "metatableAms.txt" located in the folder "Inputfiles".
# The scripts reflect the format and namegiving of the raw PWS and radar observations that were used in our example. Be aware that it needs to be rewritten to be applicable for raw time series that were stored in a different manner. This script can then serve as an example to standardize the raw data.  

## step 1: set parameters and generate neighbourlist
# Change parameter settings in InputFiles/Filtersettings.txt. Give this set of parameter settings a unique name, e.g. "A".
# Run PWSQC_1.R
# notes: the unique parameter settings label ensures that results are not overwritten. The parameter settings do not affect raw data ('Ndataset' and 'Rdataset'). Do not make changes in InputFiles/Filtersettings.txt before all steps have been completed. Any notes/comments in InputFiles/Filtersettings.txt other than changes in the parameter values should be made with precurser '#' to avoid errors. The method to calculate distances between stations is appropriate for locations in the Netherlands, but needs to be adjusted for study areas in other parts of the world due to the projection.  


  # # QC # #

## step 2: attribute FZ and HI flags
# run PWSQC_2.R
# notes: makes use of 'Ndataset' and 'neighbourlist'

## step 3: attribute SO and biascorrectionfactor
# run PWSQC_3.R
# notes: makes use of 'Ndataset', 'neighbourlist', 'FZ_flags', 'HI_flags'. As this process is very timeconsuming this script can best be started at the end of day so it can run overnight / over the weekend. 


  # # Validation # # 

## Step 4: check if the attributed flags accurately identify FZ, HI, SO and bias errors in the dataset
# run PWSQC_4.R
# notes: the script makes a graph for visual inspection of the results, as wel as calculates and prints validation metrics using 'Rdataset' as ground-truth. In the example a period of 13 months is run where we only validate the last 12 months in order to allow for a lead-up time of the filters. 



