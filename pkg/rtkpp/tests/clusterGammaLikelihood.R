# TODO: Add comment
#
# Author: iovleff
#-----------------------------------------------------------------------
library(rtkpp)
data(iris)

gamma_model <- clusterGamma(iris[1:4], nbCluster = 2:6, modelNames = gammaNames(shapeBetweenCluster = "all"), strategy = clusterStrategy(nbTry = 3, nbInit = 5))

data<-gamma_model@data
nbSample <- nrow(data)
nbVariable <- ncol(data)
nbCluster <- gamma_model@nbCluster
prop <- gamma_model@pk
shape <- gamma_model@shape
scale <- sqrt(gamma_model@scale)

sum<-0
for (i in 1:nbSample)
{
  max <- -Inf
  for (k in 1:nbCluster)
  {
    sum2 <-log(prop[k]);
    for (j in 1:nbVariable)
    {  sum2 = sum2 + dgamma(data[i,j], shape=shape[k,j], scale=scale[k,j],log=TRUE)}
    max = max(max, sum2);
  }
  sum1 <- 0
  for (k in 1:nbCluster)
  {
    sum2 <-log(prop[k]);
    for (j in 1:nbVariable)
    {  sum2 = sum2 + dgamma(data[i,j], shape=shape[k,j], scale=scale[k,j],log=TRUE)}
    sum1 = sum1 + exp(sum2-max);
  }
  sum = sum + log(sum1) + max
}
cat("Computed log-likelihood: ", sum, "\n")
cat("Model log-likelihood: ", gamma_model@lnLikelihood, "\n")
