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

dir <- integer(1) # direction of extreme (warm (1) or cold (-1))
p <- integer(1) # dimensionality of trimmed covariance matrix

calcMDsqFromMVG <- function(rectGrid, origYData, day, data)
{
    # calculates the parameters of the baseline multivariate normal that
    # corresponds to the rectangle under investigation, then calculates the 
    # mahalanobis distance of current observation from that MVG
    #
    # Args
    # rectGrid: data frame containing information about the rectangle under
    # investigation
    # origYdata: original data for all years in order to facilitate calculation
    # of covariance matrix of baseline distribution
    # day: current day under consideration
    # data: the data corresponding to the current day
    #
    # Returns
    # MDsq: Squared Mahalanobis distance of rectangle observation from MVG
    # dir: not returned, but explicitly sets direction of extreme
    
    # current observed vector under the window
    currObs <- rectGrid[, 3]
    # gather data under spatio-temporal cuboid specified by window in d
    d <- matrix(0, nrow = 175)
    for(g in 1:nrow(rectGrid))
    {
        lat <- rectGrid[g, 1]
        long <- rectGrid[g, 2]
        timeWindow <- (day - 2):(day + 2)
        y <- unlist(lapply(origYData, "[", timeWindow, lat, long)) # 175 values
        d <- cbind(d, y)
    }
    d <- d[, -1]
    # remove from d the observation that's currently under inspection for
    # anomalous behavior so that estimates aren't biased
    if(nrow(rectGrid) == 1)
    {
        d <- d[-which(d == currObs)]
        mu <- mean(d)
    } else {
        matchidx <- apply(d, 1, FUN = function(obs) all(obs == currObs))
        d <- d[-which(matchidx), ]
        mu <- colMeans(d)
    }
    # if window spans warm and cold extremes (Quadrants II or IV), skip
    if(!(all((rectGrid[, 3] - mu) >= 0) |
         all((rectGrid[, 3] - mu) <= 0)))
    {
        p <<- 0
        return(0)
    }
    # find direction of anomalous behavior
    dir <<- ifelse(mean(rectGrid[, 3]) > mean(mu), 1, -1)
    COV <- cov(as.matrix(d)) # covariance matrix
    iCOV <- pseudoinverse(COV)
    # squared mahalanobis distance
    MDsq <- mahalanobis(currObs, center = mu, cov = iCOV, inverted = TRUE) 
    return(MDsq)
}

# load the data if not already done so
if(!exists("origYData"))
    source("loadData.R")

# ten random days in the period of record
years <- c(1979)
days <- c(3)

plotpath <- paste("~/Documents/Scanning window/RandomDaysPlots/plot%02d.jpg")
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
                            rectGrid[nrow(rectGrid) + 1, ] <- c(i, j,
                                                                data[i, j])
                        }
                    }
                    p <- nrow(rectGrid) # no. of grid cells in window / dimensionality
                    MDsq <- calcMDsqFromMVG(rectGrid, origYData, day,
                                            data)
                    Pr <- ifelse(p == 0, 0, pchisq(MDsq, p))
                    # color the grid boxes with value = +-Pr
                    for(g in 1:nrow(rectGrid))
                    {
                        val <- max(abs(resToday[rectGrid[g, 1],
                                                rectGrid[g, 2]]), Pr)
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