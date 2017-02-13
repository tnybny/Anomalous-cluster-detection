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

calcMDsqFromMVG <- function(rectGrid, origYData, day)
{
    # calculates the parameters of the baseline multivariate normal that
    # corresponds to the rectangle under investigation, then calculates the 
    # mahalanobis distance of current observation from that MVG
    #
    # Params
    # rectGrid: data frame containing information about the rectangle under
    # investigation
    # origYdata: original data for all years in order to facilitate calculation
    # of covariance matrix of baseline distribution
    # day: current day under consideration
    #
    # Returns
    # MDsq: Squared Mahalanobis distance of rectangle observation from MVG
    # dir: not returned, but explicitly sets direction of extreme
    
    # observed temperatures under current window
    currObs <- rectGrid[, 3]
    # gather data under spatio-temporal cuboid specified by window in d
    timeWindow <- (day - 2):(day + 2)
    d <- t(do.call(rbind, lapply(1:nrow(rectGrid), FUN = function(g) {
        unlist(lapply(origYData, "[", timeWindow, rectGrid[g, 1],
                      rectGrid[g, 2]))})))
    # remove from d the observation that's currently under inspection for
    # anomalous behavior so that mean and covariance estimates aren't biased
    if(nrow(rectGrid) == 1)
    {
        d <- d[-which(d == currObs)[1]]
        mu <- mean(d)
    } else {
        matchidx <- apply(d, 1, FUN = function(obs) all(obs == currObs))
        d <- d[-which(matchidx)[1], ]
        mu <- colMeans(d)
    }
    # if window spans warm and cold extremes (Quadrants II or IV), skip
    if(!(all((currObs - mu) >= 0) |
         all((currObs - mu) <= 0)))
    {
        return(0)
    }
    # find direction of anomalous behavior
    dir <<- ifelse(mean(currObs) > mean(mu), 1, -1)
    iCOV <- pseudoinverse(cov(as.matrix(d)))
    # squared mahalanobis distance
    MDsq <- mahalanobis(currObs, center = mu, cov = iCOV, inverted = TRUE) 
    return(MDsq)
}

# load the data if not already done so
if(!exists("origYData"))
    source("loadData.R")

# period of record
years <- 1979
days <- 3

plotpath <- paste("./allplots/plot%02d.jpg")
jpeg(plotpath, width = 1024, height = 680)

ptm <- proc.time()

for(year in years){
    for(day in days){
        mat <- origYData[[year - 1978]][day, , ]
        
        # create result matrix to color
        resToday <- matrix (0, 73, 144)
        
        for (i1 in 1:73){ 
            for (j1 in 1:144){
                for (i2 in i1:(i1 + 2)){
                    if(i2 > 73)
                        i2 <- 73
                    for (j2 in j1:(j1 + 2)){
                        if(j2 > 144)
                            j2 <- 144
                        rectGrid <- data.frame(expand.grid(i1:i2, j1:j2),
                                               c(mat[i1:i2, j1:j2]))
                        MDsq <- calcMDsqFromMVG(rectGrid, origYData, day)
                        Pr <- ifelse(MDsq == 0, 0, pchisq(MDsq, nrow(rectGrid)))
                        # color the grid cells with value = +-Pr
                        changeIdx <- !(abs(c(resToday[i1:i2, j1:j2])) > Pr)
                        resToday[i1:i2, j1:j2][changeIdx] <- dir * Pr
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
}

print(proc.time() - ptm)

dev.off()