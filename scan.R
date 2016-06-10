# Author: Bharathkumar Ramachandra/tnybny
# this file performs the scanning window based approach to detecting
# clusters of anomalous behavior

# clear workspace
rm(list = ls(all = T))

# load required libraries
require(R.matlab) # for mat files
require(rworldmap) # for mapping
require(colorspace) # for legend
require(corpcor) # for pseudoinverse

# set the color palette
palette(diverge_hsv(21))

# load the baseline
baselinepath <- paste("~/Documents/Scanning window/Baseline/", sep = "")
meanBase <- readMat(paste(baselinepath, "meanBaseline.mat", sep = ""))$meanBase
dir <- integer(1) # direction of extreme (warm (1) or cold (-1))
p <- integer(1) # dimensionality of trimmed covariance matrix

calcMDsqFromMVG <- function(rectGrid, meanBase, origYData, day, data)
{
    # calculates the parameters of the baseline multivariate normal that
    # corresponds to the rectangle under investigation, then calculates the 
    # mahalanobis distance of current observation from that
    #
    # Args
    # rectGrid: data frame containing information about the rectangle under
    # investigation
    # meanBase: pre-computed long term mean Baseline
    # origYdata: original data for all years in order to facilitate calculation
    # of covariance matrix of baseline distribution
    # day: current day under consideration
    # data: the data corresponding to the current day
    #
    # Returns
    # MD: Mahalanobis distance of rectangle observation from baseline
    # dir: not returned, but explicitly sets direction of extreme
    x <- sapply(1:nrow(rectGrid), FUN = function(i, a, b) b[a[i, 1], a[i, 2]],
                rectGrid, data)
    mu <- sapply(1:nrow(rectGrid),
                 FUN = function(i, a, b, c) b[c, a[i, 1], a[i, 2]],
                 rectGrid, meanBase, day)
    # if window spans warm and cold extremes (Quadrants II or IV), skip
    if(!(all((rectGrid[, 3] - mu) >= 0) |
         all((rectGrid[, 3] - mu) <= 0)))
    {
        p <<- 0
        return(0)
    }
    # find direction of anomalous behavior
    dir <<- ifelse(mean(rectGrid[, 3]) > mean(mu), 1, -1)
    COV <- covBaseRectScore(rectGrid, origYData, day)
    iCOV <- pseudoinverse(COV)
    MDsq <- mahalanobis(x, center = mu, cov = iCOV, inverted = TRUE) 
    # squared mahalanobis distance
    return(MDsq)
}

covBaseRectScore <- function(rectGrid, origYData, day)
{
    # computes the covariance matrix for the baseline corresponding to the
    # rectangle under consideration
    #
    # Args
    # rectGrid: data frame containing information about the rectangle under
    # investigation
    # origYdata: original data for all years in order to facilitate calculation
    # of covariance matrix of baseline distribution
    # day: current day under consideration
    #
    # Returns
    # COV: covariance matrix
    d <- matrix(0, nrow = 175)
    for(g in 1:nrow(rectGrid))
    {
        lat <- rectGrid[g, 1]
        long <- rectGrid[g, 2]
        timeWindow <- (day - 2):(day + 2)
        x <- unlist(lapply(origYData, "[", timeWindow, lat, long)) # 175 values
        d <- cbind(d, x)
    }
    d <- d[, -1]
    COV <- cov(as.matrix(d))
    return(COV)
}

# load the data if not already done so
if(!exists("origYData"))
    source("loadData.R")

# ten random days in the period of record
years <- c(1979)
days <- c(3)

plotpath <- paste("~/Documents/Scanning window/TenRandomDaysPlots/plot%02d.jpg")
jpeg(plotpath, width = 1024, height = 680)

for(it in 1:length(days)){
    year <- years[it]
    day <- days[it]
    data <- origYData[[year - 1978]][day, , ]
    
    # create object to color
    resToday <- matrix (0, 73, 144)
    
    for (i1 in 1:73){ 
        print(i1)
        for (j1 in 1:144){
            for (i2 in i1:(i1 + 2)){
                if(i2 > 73){
                    i2 <- 73
                }
                for (j2 in j1:(j1 + 2)){
                    if(j2 > 144){
                        j2 <- 144
                    }
                    rectGrid <- data.frame(lat = numeric(), 
                                           long = numeric(),
                                           val = double())
                    for(i in i1:i2)
                    {
                        for(j in j1:j2)
                        {
                            rectGrid[nrow(rectGrid) + 1, ] <- c(i, j, data[i, j])
                        }
                    }
                    p <- nrow(rectGrid) # number of grid cells in window
                    MDsq <- calcMDsqFromMVG(rectGrid, meanBase, origYData, day,
                                            data)
                    v <- ifelse(p == 0, 0, pchisq(MDsq, p))
                    # color the grid boxes with value = +-MD
                    for(g in 1:nrow(rectGrid))
                    {
                        val <- max(abs(resToday[rectGrid[g, 1], rectGrid[g, 2]]),
                                  v)
                        resToday[rectGrid[g, 1], rectGrid[g, 2]] <- dir * val
                    }
                }
            }
        }
    }
    
    # do mapping transformations
    resToday <- t(resToday)
    resToday <- resToday[, ncol(resToday):1]
    
    # cut at 5% at both tails
    cutoff <- sort(abs(resToday), decreasing = T)[(0.05 * length(resToday))]
    cutResToday <- resToday
    cutResToday[cutResToday < cutoff & cutResToday > -cutoff] = 0
    
    # linearly stretch the values to [-1, 1] scale
    negRange <- range(cutResToday[cutResToday < 0])
    posRange <- range(cutResToday[cutResToday > 0])
    cutResTodayPrime <- ifelse(cutResToday < 0,
                               (cutResToday - negRange[1]) / 
                                   (negRange[2] - negRange[1]) *
                                   (-0.001 - (-1)) + (-1),
                               cutResToday)
    cutResTodayPrime <- ifelse(cutResToday > 0,
                               (cutResToday - posRange[1]) / 
                                   (posRange[2] - posRange[1]) *
                                   (1 - (0.001)) + (0.001),
                               cutResTodayPrime)
    
    # map the result
    mapGriddedData(cutResTodayPrime, numCats = 21, catMethod = "diverging",
                   colourPalette = "palette", borderCol = "black")
}

dev.off()