# Author: Bharathkumar Ramachandra/tnybny
# this file tests random windows for normality

# clear workspace
rm(list = ls(all = T))

# load required libraries
require(MVN)

# load the data if not already done so
if(!exists("origYData"))
    source("loadData.R")

n <- 1000 # number of windows to test
# randomly sample time and space indices
years <- sample(1979:2013, n, replace = T)
days <- sample(3:363, n, replace = T)
lats <- sample(15:58, n, replace = T)
longs <- sample(3:141, n, replace = T)

ret <- vector()
p.values.skewness <- vector()
p.values.kurtosis <- vector()
for(i in 1:n)
{
    year <- years[i]
    day <- days[i]
    lat <- lats[i]
    long <- longs[i]
    data <- origYData[[year - 1978]][day, , ]
    i1 <- lat
    j1 <- long
    i2 <- i1 + 1
    j2 <- j1 + 1
    # rectGrid is the current window under observation
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
    d <- matrix(0, nrow = 175) # d holds the data that needs to be tested for normality
    for(g in 1:nrow(rectGrid))
    {
        lat <- rectGrid[g, 1]
        long <- rectGrid[g, 2]
        timeWindow <- (day - 2):(day + 2)
        x <- unlist(lapply(origYData, "[", timeWindow, lat, long)) # 175 values
        d <- cbind(d, x)
    }
    d <- d[, -1]
    toremove <- vector()
    for(j in 1:3)
    {
        if(j %in% toremove)
            next
        for(k in (j + 1):4)
        {
            if(k %in% toremove)
                next
            if(cov(d[, j], d[, k]) > 0.8)
            {
                toremove <- c(toremove, k)
            }
        }
    }
    if(length(toremove))
        d <- d[, -toremove]
    if(length(toremove) == 3)
    {
        a <- shapiro.test(d)
        if(a$p.value <= 0.01)
            ret <- c(ret, 0) # not normal
        else
            ret <- c(ret, 1) # normal
    } else {
        a <- mardiaTest(d)
        p.values.kurtosis <- c(p.values.kurtosis, a@p.value.kurt)
        p.values.skewness <- c(p.values.skewness, a@p.value.skew)
        if(a@p.value.skew > 0.01 && a@p.value.kurt > 0.01)
            ret <- c(ret, 1) # normal
        else
            ret <- c(ret, 0) # not normal
    }
}

print(sum(ret == 1))
# passes around 55% of times for dimensionality 4 at 1% significance level