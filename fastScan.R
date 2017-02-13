# Author: Bharathkumar Ramachandra/tnybny
# this file performs the scanning window based approach to detecting
# clusters of anomalous behavior - sacrifice space for time.
# exactly the same as scan.R but faster, except for 10+ decimal differences.

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

# load the data if not already done so
if(!exists("origYData"))
    source("loadData.R")

# create result structure to color
res <- list()
for(i in 1:35)
    res[[i]] <- array(0, c(365, 73, 144))

ptm <- proc.time()

years <- 1:35
days <- 3

for (i1 in 1:73){
    print(i1)
    for (j1 in 1:144){
        for (i2 in i1:(i1 + 2)){
            if(i2 > 73)
                i2 <- 73
            for (j2 in j1:(j1 + 2)){
                if(j2 > 144)
                    j2 <- 144
                rectGrid <- data.frame(expand.grid(i1:i2, j1:j2))
                for(day in days)
                {
                    timeWindow <- (day - 2):(day + 2)
                    d <- t(do.call(rbind, lapply(1:nrow(rectGrid),
                                                 FUN = function(g) {
                        unlist(lapply(origYData, "[", timeWindow,
                                      rectGrid[g, 1], rectGrid[g, 2]))})))
                    for(year in years)
                    {
                        # remove from d the observation that's currently under inspection for
                        # anomalous behavior so that mean and covariance estimates aren't biased
                        currObs <- c(origYData[[year]][day, i1:i2, j1:j2])
                        if(nrow(rectGrid) == 1)
                        {
                            d <- d[-which(d == currObs)[1]]
                            mu <- mean(d)
                        } else {
                            matchidx <- apply(d, 1, FUN = function(obs) 
                                all(obs == currObs))
                            d <- d[-which(matchidx)[1], ]
                            mu <- colMeans(d)
                        }
                        # if window spans warm and cold extremes 
                        # (Quadrants II or IV), skip
                        if(!(all((currObs - mu) >= 0) |
                             all((currObs - mu) <= 0)))
                            next
                        # find direction of anomalous behavior
                        dir <- ifelse(mean(currObs) > mean(mu), 1, -1)
                        iCOV <- pseudoinverse(cov(as.matrix(d)))
                        # squared mahalanobis distance
                        MDsq <- mahalanobis(currObs, center = mu, cov = iCOV,
                                            inverted = TRUE)
                        Pr <- ifelse(MDsq == 0, 0, pchisq(MDsq, nrow(rectGrid)))
                        # color the grid cells with value = +-Pr
                        changeIdx <- !(abs(c(res[[year]][day, i1:i2, j1:j2])) >
                                           Pr)
                        res[[year]][day, i1:i2, j1:j2][changeIdx] <- dir * Pr
                    }
                }
            }
        }
    }
}

print(proc.time() - ptm)

plotpath <- paste("./fastplots/plot%02d.jpg")
jpeg(plotpath, width = 1024, height = 680)

for(year in years)
{
    for(day in days)
    {
        resToday <- t(res[[year]][day, , ])
        resToday <- resToday[, ncol(resToday):1]
        
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

dev.off()