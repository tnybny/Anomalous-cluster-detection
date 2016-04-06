# Author: Bharathkumar Ramachandra/tnybny
# this file performs the scanning window based approach to detecting
# clusters of anomalous behavior

# clear workspace
rm(list = ls(all = T))

# load required libraries
require(R.matlab)
require(rworldmap)
require(colorspace)
# set the color palette
palette(diverge_hsv(21))

# load the baseline
baselinepath <- "./CSVs/Scanning_window/Baseline/"
meanBase <- readMat(paste(baselinepath, "meanBaseline.mat", sep = ""))$meanBase
MDthresh <- double(1)
dir <- double(1)
toRemove <- numeric()
p <- integer(1)

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
    COV = covBaseRectScore(rectGrid, origYData, day)
    mu <- double()
    x <- double()
    for(i in 1:nrow(rectGrid))
    {
        x <- c(x, data[rectGrid[i, 1], rectGrid[i, 2]])
        mu <- c(mu, meanBase[day, rectGrid[i, 1], rectGrid[i, 2]])
    }
    if(length(toRemove) != 0)
    {
        x <- x[-toRemove]
        mu <- mu[-toRemove]
    }
    p <<- length(mu)
    # if window spans warm and cold extremes (Quadrants II or IV), skip
    if(!(all((rectGrid[-toRemove, 3] - mu) >= 0) |
             all((rectGrid[-toRemove, 3] - mu) <= 0)))
    {
        p <<- 0
        return(0)
    }
    # find direction of anomalous behavior
    if(mean(rectGrid[, 3]) > mean(mu))
        dir <<- 1
    else
        dir <<- -1
    MDsq <- mahalanobis(x, center = mu, cov = COV) # squared mahalanobis distance
#     if(length(mu) > 1)
#     {
#         MDsq <- MDsq / dim(COV)[2]   # correct for dimensionality effect
#     }
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
    # toRemove: not returned, but sets the grid boxes to remove from 
    # consideration because of multicollinearity
    d <- matrix(, nrow = 35)
    for(g in 1:nrow(rectGrid))
    {
        lat = rectGrid[g, 1]
        long = rectGrid[g, 2]
        x <- double()
        for(year in 1:35)
        {
            x <- c(x, origYData[[year]][day, lat, long])
        }
        d <- cbind(d, x)
    }
    d <- d[, -1]
    if(is.vector(d))
    {
        COV <- var(d)
    } else {
        toRemove <<- numeric()
        if(det(cov(d)) < 10 ^ (-5))
        {
            # account for perfect multicollinearity
            for(k in 1:(ncol(d) - 1))
            {
                for(l in (k + 1):ncol(d))
                {
                    if(cor(d[, k], d[, l]) > 0.95 & !(l %in% toRemove))
                    {
                        toRemove <<- c(toRemove, l)
                    }
                }
            }
            d <- d[, -toRemove]
        }
        if(is.atomic(d) || is.list(d))
        {
            COV <- var(d)
        } else {
            COV <- cov(d)
        }
    }
    return(COV)
}

# read in all the STemp data
origpath <- "./Original_Stemp_and_Z500_data/STemp_original_data/R\ data/"
origYData <- list()
for(i in 1:35)
{
    print(i)
    filename <- paste(origpath, "ncep1.STemp.", 1978 + i, ".R.mat", sep = "")
    origYData[[i]] <- readMat(filename)$STemp.dy
    if(dim(origYData[[i]])[1] == 366)
    {
        origYData[[i]] <- origYData[[i]][-366, , ]
    }
}

# for each year
year = 1979

# for each day
day = 3
data <- origYData[[year - 1978]][day, , ]
dataM1 <- origYData[[year - 1978]][day - 1, , ]  # minus 1
dataM2 <- origYData[[year - 1978]][day - 2, , ]
dataP1 <- origYData[[year - 1978]][day + 1, , ]
dataP2 <- origYData[[year - 1978]][day + 2, , ]  # plus 2

# create object to color
resToday <- matrix (0, 73, 144)

for (i1 in 1:73){ 
    print(i1)
    for (j1 in 1:144){
        for (i2 in i1:(i1 + 2)){
            if(i2 > 73){
                i2 = 73
            }
            for (j2 in j1:(j1 + 2)){
                if(j2 > 144){
                    j2 = 144
                }
                rectGrid <- data.frame(lat = numeric(), long = numeric(),
                                       val = double())
                for(i in i1:i2)
                {
                    for(j in j1:j2)
                    {
                        rectGrid[nrow(rectGrid) + 1, ] = c(i, j, data[i, j])
                        rectGrid[nrow(rectGrid) + 1, ] = c(i, j, dataM1[i, j])
                        rectGrid[nrow(rectGrid) + 1, ] = c(i, j, dataM2[i, j])
                        rectGrid[nrow(rectGrid) + 1, ] = c(i, j, dataP1[i, j])
                        rectGrid[nrow(rectGrid) + 1, ] = c(i, j, dataP2[i, j])
                    }
                }
                MDsq <- calcMDsqFromMVG(rectGrid, meanBase, origYData, day,
                                      data)
                v = ifelse(p == 0, 0, pchisq(MDsq, p))
                # color the grid boxes with value = +-MD
                for(g in 1:nrow(rectGrid))
                {
                    val = max(abs(resToday[rectGrid[g, 1], rectGrid[g, 2]]),
                              v)
                    resToday[rectGrid[g, 1], rectGrid[g, 2]] = dir * val
                }
            }
        }
    }
}

# do mapping transformations
resToday <- t(resToday)
resToday <- resToday[, ncol(resToday):1]
resToday <- resToday * 100
resToday <- round(resToday, 1)

# map the result
mapGriddedData(resToday, numCats = 21, catMethod = "diverging",
               colourPalette = "palette", borderCol = "black")

# cut at 5% at both tails
cutoff <- sort(abs(resToday), decreasing = T)[(0.05 * length(resToday))]
cutResToday <- resToday
cutResToday[cutResToday < cutoff & cutResToday > -cutoff] = 0

# map the result
mapGriddedData(cutResToday, numCats = 21, catMethod = "diverging",
               colourPalette = "palette", borderCol = "black")

# write csv
#write.table(resToday, filename, row.names = F, col.names = F, sep = ",")

# > mahalanobis(c(2.2567, 2.2567), center = c(0, 0), cov = matrix(c(1, 0.7, 0.7, 1), 2, 2))
# [1] 5.991406
# > qchisq(0.95, 2)
# [1] 5.991465
# > pchisq(5.991465, 2) # MD^2 is chi-squared distributed!
# [1] 0.95