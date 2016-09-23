# clear workspace
rm(list = ls (all = T))

# load required librarues
library(ggplot2)
library(MASS)

# create 2 1-d rnorms with no correlation
x1 <- rnorm(2000, 0, 1)
x2 <- rnorm(2000, 0, 1)
d <- data.frame(x1, x2)

# plot the points on each axis separately
g <- ggplot() + theme(text = element_text(size=20)) +
    geom_point(data = d, aes(x = x1, y = 0), col = "red", size = 0.5) +
    geom_point(data = d, aes(x = 0, y = x2), col = "blue", size = 0.5) 

# add 95% intervals for each axis
g <- g + geom_segment(aes(x = -1.96, y = -0.1, xend = -1.96, yend = 0.1)) + 
    geom_segment(aes(x = 1.96, y = -0.1, xend = 1.96, yend = 0.1)) +
    geom_segment(aes(x = -0.1, y = -1.96, xend = 0.1, yend = -1.96)) + 
    geom_segment(aes(x = -0.1, y = 1.96, xend = 0.1, yend = 1.96)) + ylab("x2")

# add 95% confidence ellipse andn title
g <- g + stat_ellipse(type = "norm", data = d, aes(x = x1, y = x2)) +
    coord_fixed() 
g

# now add some correlation
m <- data.frame(mvrnorm(n = 2000, mu = c(0, 0),
                        Sigma = matrix(c(1, 0.7, 0.7, 1), 2, 2)))
names(m) <- c("x1", "x2")

# plot the points in 2 dimensions
h <- ggplot(data = m, aes(x = x1, y = x2)) + 
    theme(text = element_text(size=20)) +
    geom_point(size = 0.5, col = "orange")

# add 95% intervals for each axis
h <- h + geom_segment(aes(x = -2, y = -0.1, xend = -2, yend = 0.1), col = "red") + 
    geom_segment(aes(x = 2, y = -0.1, xend = 2, yend = 0.1), col = "red") +
    geom_segment(aes(x = -0.1, y = -2, xend = 0.1, yend = -2), col = "blue") + 
    geom_segment(aes(x = -0.1, y = 2, xend = 0.1, yend = 2), col = "blue") +
    ylab("x2") + coord_fixed()

# add 95% confidence ellipse and title
h <- h + stat_ellipse(type = "norm") +
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0)

# add demo points
h <- h + 
    geom_point(x = 2.08, y = 1.9, color = "green", size = 4, shape = 18) + 
    geom_point(x = -0.1, y = -1.90, color = "magenta", size = 4, shape = 18)
h + annotate("text", x = 1.8, y = 1.9, label = "1", col = "red", size = 8) +
    annotate("text", x = -0.4, y = -1.8, label = "2", col = "red", size = 8)
