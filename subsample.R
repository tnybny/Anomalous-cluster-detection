#resample

rescale <- function(x, newrange = range(x)){
    xrange <- range(x)
    mfac <- (newrange[2] - newrange[1]) / (xrange[2] - xrange[1])
    newrange[1] + (x - xrange[1]) * mfac
}

ResizeMat <- function(mat, ndim=dim(mat)){
    if(!require(fields)) stop("`fields` required.")
    
    # input object
    odim <- dim(mat)
    obj <- list(x = 1:odim[1], y = 1:odim[2], z = mat)
    
    # output object
    ans <- matrix(NA, nrow = ndim[1], ncol = ndim[2])
    ndim <- dim(ans)
    
    # rescaling
    ncord <- as.matrix(expand.grid(seq_len(ndim[1]), seq_len(ndim[2])))
    loc <- ncord
    loc[, 1] = rescale(ncord[, 1], c(1, odim[1]))
    loc[, 2] = rescale(ncord[, 2], c(1, odim[2]))
    
    # interpolation
    ans[ncord] <- interp.surface(obj, loc)
    
    ans
}