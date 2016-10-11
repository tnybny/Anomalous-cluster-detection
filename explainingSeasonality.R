# Author: Bharathkumar Ramachandra/tnybny

# clear workspace
rm(list = ls(all = T))

# load required libraries
require(MVN)

# load the data if not already done so
if(!exists("origYData"))
    source("loadData.R")

lat <- 55
long <- 126

plot(x = NULL, y = NULL, xlim = c(1, 365),
     ylim = c(min(origYData[[1]][, lat, long]) - 2,
              max(origYData[[1]][, lat, long]) + 2),
     type = 'n', xlab = "Day of year", ylab = "Near surface temperature",
     mgp = c(1.5, 0.6, 0))
for(i in 1:5)
{
    data <- origYData[[i]][, lat, long]
    lines(data, col = i)
}
lines(c(100, 100, 105, 105, 100), c(262, 290, 290, 262, 262),
      lwd = 2, col = "black")

plot(x = NULL, y = NULL, xlim = c(1, 5),
     ylim = c(min(origYData[[1]][(102 - 2):(102 + 2), lat, long]) - 5,
              max(origYData[[1]][(102 - 2):(102 + 2), lat, long]) + 10),
     type = 'n', xlab = "Day in temporal window", ylab = "Temperature")
data <- do.call(cbind,
                lapply(1:35, function(i) origYData[[i]][(102 - 2):(102 + 2),
                                                         lat, long]))
means <- rowMeans(data)
sds <- apply(data, 1, sd)
lines(means)
lines(means + 2 * sds, col = "red")
lines(means - 2 * sds, col = "red")

data <- as.vector(data)
hist(data, prob = TRUE, xlab = "Temperature", main = "")
lines(density(data))
a <- shapiro.test(data)
text(285, 0.1,
     labels = paste("p value from Shapiro-Wilk test\n = ", round(a$p.value, 4),
                    sep = ""))