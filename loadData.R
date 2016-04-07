# Author: Bharathkumar Ramachandra/tnybny
# this file jsut reads in all the temperature data

# read in all the STemp data
origpath <- paste("~/Google Drive/NCSU/STAC lab research/",
                  "Monster Ridges project/Data/Original_Stemp_and_Z500_data/",
                  "STemp_original_data/R\ data/", sep = "")
origYData <- list()
for(i in 1:35)
{
    filename <- paste(origpath, "ncep1.STemp.", 1978 + i, ".R.mat", sep = "")
    origYData[[i]] <- readMat(filename)$STemp.dy
    if(dim(origYData[[i]])[1] == 366)
    {
        origYData[[i]] <- origYData[[i]][-366, , ] # just drop the last day for
        # leap years
    }
}