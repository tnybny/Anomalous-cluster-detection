# scan on natural images
rm(list = ls(all = T))

require(tiff) # 
require(colorspace) # for legend
require(corpcor) # for pseudoinverse
require(raster) # for downsampling matrix

dirs <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE)
validDirs <- dirs[nchar(dirs) == 9]
files <- lapply(validDirs, list.files, full.names = TRUE)
rm(dirs, validDirs)

source("~/Documents/Probabilistic Anomalous Cluster Detection/subsample.R")
imgs <- list()

for(i in 1:16)
    imgs[[i]] <- simplify2array(lapply(files[[i]], readTIFF))

# uncomment below to subsample the matrix to half size
# for(i in 1:16)
#     imgs[[i]] <- simplify2array(lapply(lapply(files[[i]], readTIFF),
#                                        ResizeMat, c(79, 119)))

jpeg("Image%03d.jpg", quality = 100)

for(v in 5)
{
    for(t in 23)
    {
        print(t)
        res <- array(0, dim = dim(imgs[[1]][, , t]))
        mat <- imgs[[v]][, , t]
        for (i1 in 1:dim(mat)[1]){
            print(i1)
            for (j1 in 1:dim(mat)[2]){
                for (i2 in i1:(i1 + 5)){
                    if(i2 > dim(mat)[1])
                        i2 <- dim(mat)[1]
                    for (j2 in j1:(j1 + 5)){
                        if(j2 > dim(mat)[2])
                            j2 <- dim(mat)[2]
                        if((i2 - i1) * (j2 - j1) < 2 | (i2 - i1) * (j2 - j1) > 16)
                            next
                        rectGrid <- data.frame(expand.grid(i1:i2, j1:j2),
                                               c(mat[i1:i2, j1:j2]))
                        p <- nrow(rectGrid)
                        d <- sapply(1:nrow(rectGrid),
                               function(x) {imgs[[v]][rectGrid[x, 1],
                                                      rectGrid[x, 2], 1:200]})
                        currObs <- rectGrid[, 3]
                        matchidx <- apply(d, 1,
                                          FUN = function(obs) all(
                                              obs == currObs))
                        d <- d[-matchidx, ]
                        mu <- colMeans(d)
                        if(!(all((currObs - mu) >= 0) |
                             all((currObs - mu) <= 0)))
                            next
                        iCOV <- pseudoinverse(cov(as.matrix(d)))
                        MDsq <- mahalanobis(currObs, center = mu, cov = iCOV,
                                            inverted = TRUE) 
                        Pr <- ifelse(p == 0, 0, pchisq(MDsq, nrow(rectGrid)))
                        # color the grid cells with value = Pr
                        changeIdx <- !(abs(c(res[i1:i2, j1:j2])) > Pr)
                        res[i1:i2, j1:j2][changeIdx] <- Pr
                    }
                }
            }
        }
        res <- t(apply(res, 2, rev))
        image(res, col = grey(seq(0, 1, length = 256)))
        
        # uncomment below for 2% anomaly cutoff
        # cutoff <- sort(res, decreasing = T)[(0.02 * length(res))]
        # cutRes <- res
        # cutRes[cutRes < cutoff & cutRes > -cutoff] = 0
        # 
        # image(cutRes, col = grey(seq(0, 1, length = 256)))
    }
}

dev.off()