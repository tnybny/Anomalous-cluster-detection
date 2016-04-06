# this is to investigate the claim that sum of squared MDs
# equals rank(covariance matrix) * (n - 1) where n is the number of points used
# to generate the MVG

library(MASS)
library(Matrix)

m <- data.frame(mvrnorm(n = 200, mu = c(0, 0, 0),
                        Sigma = matrix(c(1, 0.7, 0.3, 0.7, 1, 0.4, 0.3, 0.4, 1),
                                       3, 3, byrow = TRUE)))

mn = apply(m, 2, mean)
sigma = cov(m)

mdsq <- vector()
for(i in 1:nrow(m))
{
    mdsq[i] <- mahalanobis(m[i, ], center = mn, cov = sigma)
}

sum(mdsq) == rankMatrix(cov(m))[1] * (nrow(m) - 1)
# its true!!


Posdef <- function (n, ev = runif(n, 0, 10)) 
{
    # function to randomly generate pos. def. matrix
    #
    # Args:
    # n - number of linearly independent dimensions of matrix
    # ev - positive eigenvalues
    #
    # Returns:
    # Pos. def. matrix with specified eigenvalues and number of dimensions
    Z <- matrix(ncol=n, rnorm(n^2))
    decomp <- qr(Z)
    Q <- qr.Q(decomp) 
    R <- qr.R(decomp)
    d <- diag(R)
    ph <- d / abs(d)
    O <- Q %*% diag(ph)
    Z <- t(O) %*% diag(ev) %*% O
    return(Z)
}

mat <- Posdef(9)
m <- data.frame(mvrnorm(n = 2000, mu = rep(0, 9),
                        Sigma = mat))
mn = apply(m, 2, mean)
sigma = cov(m)

mdsq <- vector()
for(i in 1:nrow(m))
{
    mdsq[i] <- mahalanobis(m[i, ], center = mn, cov = sigma)
}

sum(mdsq) == rankMatrix(cov(m))[1] * (nrow(m) - 1)
