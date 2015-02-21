# TODO: Add comment
#
# Author: iovleff
#-----------------------------------------------------------------------
library(rtkpp)
data(iris)

gauss_model <- clusterDiagGaussian(iris[1:4], nbCluster = 3:4, modelNames = clusterDiagGaussianNames()
              , strategy = clusterFastStrategy(), criterion = "ICL")

data<-gauss_model@component@data
nbSample <- nrow(data)
nbVariable <- ncol(data)
nbCluster <- gauss_model@nbCluster
prop <- gauss_model@pk
mean <- gauss_model@component@mean
sigma <- gauss_model@component@sigma

f <-vector(length=nbSample)
lnComp <- vector(length=nbCluster)

for (i in 1:nbSample)
{
  for (k in 1:nbCluster)
  { lnComp[k] = log(prop[k]) + sum(dnorm(data[i,], mean[k,], sigma[k,],log=TRUE)); }
  lmax <- max(lnComp)

  lnComp =  lnComp -lmax;

  f[i] = log(sum(exp(lnComp))) + lmax;
}

cat("Computed log-likelihood: ", sum(f), "\n")
cat("Model log-likelihood: ", gauss_model@lnLikelihood, "\n")
