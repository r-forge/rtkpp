# TODO: Add comment
#
# Author: iovleff
#-----------------------------------------------------------------------
library(rtkpp)
data(iris)

gauss_model <- clusterDiagGaussian(iris[1:4], nbCluster = 2:6, modelNames = diagGaussianNames(), strategy = clusterStrategy(nbTry = 3, nbInit = 5))

data<-gauss_model@data
nbSample <- nrow(data)
nbVariable <- ncol(data)
nbCluster <- gauss_model@nbCluster
prop <- gauss_model@pk
mean <- gauss_model@mean
sigma <- sqrt(gauss_model@sigma2)

sum<-0
for (i in 1:nbSample)
{
  max <- -Inf
  for (k in 1:nbCluster)
  {
    sum2 <-log(prop[k]);
    for (j in 1:nbVariable)
    { sum2 = sum2 + dnorm(data[i,j], mean[k,j], sigma[k,j],TRUE)}
    max = max(max, sum2);
  }
  sum1 <- 0
  for (k in 1:nbCluster)
  {
    sum2 <-log(prop[k]);
    for (j in 1:nbVariable)
    {  sum2 = sum2 + dnorm(data[i,j], mean[k,j], sigma[k,j],TRUE)}
    sum1 = sum1 + exp(sum2-max);
  }
  sum = sum + log(sum1)+max
}
cat("Computed log-likelihood: ", sum, "\n")
cat("Model log-likelihood: ", gauss_model@lnLikelihood, "\n")
