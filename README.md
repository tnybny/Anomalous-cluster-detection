# Probabilistic anomalous cluster detection in gridded data

This project will focus on an approach to detect anomalous clusters in gridded data by accounting for spatial and temporal autocorrelations.

mahalanobisDistCorrection.R explores the property of sum of squared Mahalanobis distances being equal to p\*(n-1) where p is the dimensionality and n is the number of data points. This property is not essential to the procedure.

loadData.R loads all the NCEP-1 Reanalysis[1] surface temperature files into an object called origYData.

bivarate\_gaussian\_behavior.R attempts to highlight the effect of adding dimensions to a Gaussian. It also shows how a point could be non-anomalous by itself but anomalous when considered along with a spatial neighbor.

scan.R performs the scanning window based approach to anomalous cluster detection on NCEP-1 Reanalysis data [1].

testNormality.R selects a set of random points in space and time, considers a window that would be enumerated during the scanning window procedure and tests the observations that fall under it for normality using either Shapiro-Wilk test or Mardia's test. This is to prove empirically that a good number of windows satisfy the normality assumption. 

explainingSeasonality.R attemps to highlight that a 5-day temporal window is almost flat even though there is high seasonality through the year, allowing the temperature fluctuations around a given day to be modeled roughly as a Gaussian.

## Collaborators:
tnybny

bcdutton

## Contact:
tnybny@gmail.com

##References:
1) Kalnay, Eugenia, et al. "The NCEP/NCAR 40-year reanalysis project." Bulletin of the American meteorological Society 77.3 (1996): 437-471.
