# Anomalous cluster detection using a scanning window and Mahalanobis distance

This project will focus on an approach to detect anomalous clusters in gridded climate data by accounting for spatial and temporal autocorrelations.

baselineCreation.R pre-calculates 5-day long-term means and writes them to file.

mahalanobisDistCorrection.R explores the property of sum of squared Mahalanobis distances being equal to p*(n-1) where p is the dimensionality and n is the number of data points.

loadData.R loads all the NCEP-1 Reanalysis[1] surface temperature files into an object called origYData.

bivarate_gaussian_behavior.R attempts to highlight the effect of adding dimensions to a Gaussian. It also shows how a point could be non-anomalous in a lower dimension but anomalous in a higher dimensionality.

scan.R performs the scanning window based approach to anomalous cluster detection in NCEP-1 Reanalysis data [1]. 

## Collaborators:
tnybny

bcdutton

## Contact:
tnybny@gmail.com

##References:
1) Kalnay, Eugenia, et al. "The NCEP/NCAR 40-year reanalysis project." Bulletin of the American meteorological Society 77.3 (1996): 437-471.
