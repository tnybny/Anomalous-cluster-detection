## Author: Bharathkumar Ramachandra/tnybny
# this file creates the baseline for the scanning window based approach
# to detecting clusters of anomalous behavior

# clear workspace
rm(list = ls(all = T))

# load required libraries
require(R.matlab)

data <- list()

# load all the data
for(i in 1:35)
{
    print(i)
    filename <- paste("./ncep1.STemp.", 1978 + i, ".R.mat", sep = "")
    data[[i]] <- readMat(filename)$STemp.dy
    if(dim(data[[i]])[1] == 366)
    {
        data[[i]] <- data[[i]][-366, , ]
    }
}

# calculate long term means for each day, each location
meanBaseline <- apply(simplify2array(data), 1:3, mean)

# calculate 5-day window means for each day
for(d in 3:363)
{
    L <- list(meanBaseline[d - 2, , ], meanBaseline[d - 1, , ],
         meanBaseline[d, , ], meanBaseline[d + 1, , ],
         meanBaseline[d + 2, , ])
    meanBaseline[d, , ] <- apply(simplify2array(L), 1:2, mean)
}

# for days 1, 2, 364, 365
d = 1
L <- list(meanBaseline[d, , ], meanBaseline[d + 1, , ],
          meanBaseline[d + 2, , ])
meanBaseline[d, , ] <- apply(simplify2array(L), 1:2, mean)

d = 2
L <- list(meanBaseline[d - 1, , ], meanBaseline[d, , ], meanBaseline[d + 1, , ],
          meanBaseline[d + 2, , ])
meanBaseline[d, , ] <- apply(simplify2array(L), 1:2, mean)

d = 364
L <- list(meanBaseline[d - 2, , ], meanBaseline[d - 1, , ], meanBaseline[d, , ],
          meanBaseline[d + 1, , ])
meanBaseline[d, , ] <- apply(simplify2array(L), 1:2, mean)

d = 365
L <- list(meanBaseline[d - 2, , ], meanBaseline[d - 1, , ], meanBaseline[d, , ])
meanBaseline[d, , ] <- apply(simplify2array(L), 1:2, mean)

# write mean to file
writeMat(con = paste("~/Documents/Scanning window/",
                                       "Baseline/meanBaseline.mat", sep = ""),
         meanBase = meanBaseline)