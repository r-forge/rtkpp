library(MixAll)

testPredict<-function(nbTrain , nbTest)
{
  ## test categorical predictions
  train1 <- matrix( c( sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.05,0.05,0.9))
      , sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.9,0.05,0.05))
      , sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.05,0.05,0.9))
      , sample(1:3,size=nbTrain,replace=TRUE, prob = c(0.9,0.05,0.05))
    )
    , ncol =2
  )
  test1 <- matrix( c( sample(1:3,size=nbTest,replace=TRUE, prob = c(0.05,0.05,0.9))
      , sample(1:3,size=nbTest,replace=TRUE, prob = c(0.9,0.05,0.05))
      , sample(1:3,size=nbTest,replace=TRUE, prob = c(0.05,0.05,0.9))
      , sample(1:3,size=nbTest,replace=TRUE, prob = c(0.9,0.05,0.05))
    )
    , ncol =2
  )
  model <- clusterCategorical(train1,2,models = "categorical_p_pjk")
  pred  <- clusterPredict(test1,model)
  # more than 5 classification errors is abnormal
  if (abs(sum(pred@zi) - nbTest)>5) return(FALSE)
  
  ##------------------------------------------------------------------------------
  ## test Poisson predictions
  train2 <- matrix( c( rpois(nbTrain,lambda = 1), rpois(nbTrain,lambda = 10)
      , rpois(nbTrain,lambda = 1), rpois(nbTrain,lambda = 10))
    , ncol =2
  )
  test2 <- matrix( c( rpois(nbTest,lambda = 1), rpois(nbTest,lambda = 10)
      , rpois(nbTest,lambda = 1), rpois(nbTest,lambda = 10))
    , ncol =2
  )
  model <- clusterPoisson(train2,2,models = "poisson_p_lk")
  pred  <- clusterPredict(test2,model)
  # more than 5 classification errors is abnormal
  if (abs(sum(pred@zi) - nbTest)>5) return(FALSE)
  
  ##------------------------------------------------------------------------------
  ## test Gaussian predictions
  train3 <- matrix( c( rnorm(nbTrain, mean = 1, sd=1), rnorm(nbTrain,mean = 10, sd=1)
      , rnorm(nbTrain, mean = 1, sd=1), rnorm(nbTrain,mean = 10, sd=1))
    , ncol =2
  )
  test3 <- matrix( c( rnorm(nbTest,mean = 1, sd=1), rnorm(nbTest,mean = 10, sd=1)
      , rnorm(nbTest,mean = 1, sd=1), rnorm(nbTest,mean = 10, sd=1))
    , ncol =2
  )
  model <- clusterDiagGaussian(train3,2,models = "gaussian_p_s")
  pred  <- clusterPredict(test3,model)
  # more than 5 classification errors is abnormal
  if (abs(sum(pred@zi) - nbTest)>5) return(FALSE)
  
  ##------------------------------------------------------------------------------
  ## test gamma predictions
  train4 <- matrix( c( rgamma(nbTrain, shape = 1, scale=1), rgamma(nbTrain,shape = 10, scale=1)
      , rgamma(nbTrain, shape = 1, scale=1), rgamma(nbTrain,shape = 10, scale=1))
    , ncol =2
  )
  test4 <- matrix( c( rgamma(nbTest,shape = 1, scale=1), rgamma(nbTest,shape = 10, scale=1)
      , rgamma(nbTest,shape = 1, scale=1), rgamma(nbTest,shape = 10, scale=1))
    , ncol =2
  )
  model <- clusterGamma(train4, 2, models = "gamma_p_ak_b")
  pred  <- clusterPredict(test4,model)
  # more than 5 classification errors is abnormal
  if (abs(sum(pred@zi) - nbTest)>5) return(FALSE)
  
  ##------------------------------------------------------------------------------
  ## test mixed data predictions
  train <- list(train1, train2, train3, train4)
  test <- list(test1, test2, test3, test4)
  models <- c("categorical_p_pjk", "poisson_p_lk",  "gaussian_p_s","gamma_p_ak_b")
  
  model <- clusterMixedData(train, models, 2)
  pred  <- clusterPredict(test,model)
  # more than 5 classification errors is abnormal
  if (abs(sum(pred@zi) - nbTest)>5) return(FALSE)
  
  ##------------------------------------------------------------------------------
  return(TRUE)  
}

testPredict(1000, 20)
