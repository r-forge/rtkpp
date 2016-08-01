#-----------------------------------------------------------------------
#     Copyright (C) 2012-2016  Serge Iovleff, University Lille 1, Inria
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as
#    published by the Free Software Foundation; either version 2 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this program; if not, write to the
#    Free Software Foundation, Inc.,
#    59 Temple Place,
#    Suite 330,
#    Boston, MA 02111-1307
#    USA
#
#    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
#
#-----------------------------------------------------------------------
#' @include ClusterDiagGaussian.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of a learn mixture model
#'
#' This function learn the optimal mixture model when the class labels are known
#' according to the \code{criterion} among the list of model given in \code{models}.
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the estimation process.
#' @param labels vector or factors giving the label class.
#' @param models [\code{\link{vector}}] of model names to run. By default all
#' models are estimated.
#' @param prop [\code{\link{vector}}] with the proportions of each class.
#' If NULL the proportions will be estimated using the labels.
#' @param algo character defining the algo to used in order to learn the model.
#' Possible values: "simul" (default), "impute" (faster but can produce biased results).
#' @param nbIter integer giving the number of iterations to do.
#' algo is "impute" this is the maximal authorized number of iterations. Default is 100.
#' @param epsilon real giving the variation of the log-likelihood for stopping the
#' iterations. Not used if algo is "simul". Default value is 1e-08.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL" (default).
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @examples
#' 
#' ## A quantitative example with the famous iris data set
#' data(iris)
#' 
#' ## get data and target
#' x <- as.matrix(iris[,1:4]);
#' z <- as.vector(iris[,5]);
#' n <- nrow(x); p <- ncol(x);
#' 
#' ## add missing values at random
#' indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2);
#' x[indexes] <- NA;
#' 
#' ## learn model
#' model <- learnDiagGaussian( data=x, labels= z, prop = c(1/3,1/3,1/3)
#'                           , models = clusterDiagGaussianNames(prop = "equal")
#'                           )
#' 
#' ## get summary
#' summary(model)
#' 
#' ## use graphics functions
#' \dontrun{
#' plot(model)
#' }
#' ## print model
#' \dontrun{
#' print(model)
#' }
#' 
#' ## get estimated missing values
#' missingValues(model)
#'
#' @return An instance of a learned mixture model class.
#' @rdname learners
#' @aliases learnDiagGaussian
#' @author Serge Iovleff
#' @export
#'
learnDiagGaussian <- function( data, labels, prop = NULL
                             , models=clusterDiagGaussianNames(prop = "equal")
                             , algo="simul", nbIter = 100, epsilon = 1e-08
                             , criterion="ICL"
                             , nbCore = 1)
{
  # check data
  labels = .checkDataInLearner(data, labels)
  # check proportions
  prop <- .checkPropInLearner(labels, prop)
  
  # check nbIter
  if( nbIter < 1) { stop("nbIter must be greater or equal to 1")}
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?learnDiagGaussian for the list of valid criterion")}
  
  # check models
  if (!clusterValidDiagGaussianNames(models))
  { stop("models is not valid. See ?clusterDiagGaussianNames for the list of valid model names")}
  # Create model
  model = .createMixtureModel("ClusterDiagGaussian", data, labels, prop);
  # Create algorithm
  algo = learnAlgo( algo, nbIter, epsilon)
  # start estimation of the models
  resFlag = .Call("learnMixture", model, models, algo, criterion, nbCore, PACKAGE="MixAll");
  # set names
  if (resFlag != TRUE ) {cat("WARNING: An error occur during the learning process");}
  colnames(model@component@mean)  <- colnames(model@component@data);
  colnames(model@component@sigma) <- colnames(model@component@data);
  model
}

#' @rdname learners
#' @aliases learnPoisson
learnPoisson <- function( data, labels, prop = NULL
                        , models=clusterPoissonNames(prop = "equal")
                        , algo="simul", nbIter = 100, epsilon = 1e-08
                        , criterion="ICL"
                        , nbCore = 1)
{
  # check data
  labels = .checkDataInLearner(data, labels)
  # check proportions
  prop <- .checkPropInLearner(labels, prop)
  
  # check nbIter
  if( nbIter < 1) { stop("nbIter must be greater or equal to 1")}
  
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?learnPoisson for the list of valid criterion")}
  
  # check models
  if (!clusterValidPoissonNames(models))
  { stop("models is not valid. See ?clusterPoissonNames for the list of valid model names")}
  # Create model
  model = .createMixtureModel("ClusterPoisson", data, labels, prop);
  # Create algorithm
  algo = learnAlgo( algo, nbIter, epsilon)
  # start estimation of the models
  resFlag = .Call("learnMixture", model, models, algo, criterion, nbCore, PACKAGE="MixAll");
  # set names
  if (resFlag != TRUE ) {cat("WARNING: An error occur during the learning process");}
  colnames(model@component@lambda)  <- colnames(model@component@data);
  model
}

#' @rdname learners
learnGamma <- function( data, labels, prop = NULL
                      , models=clusterGammaNames(prop = "equal")
                      , algo="simul", nbIter = 100, epsilon = 1e-08
                      , criterion="ICL"
                      , nbCore = 1)
{
  # check data
  labels = .checkDataInLearner(data, labels)
  # check proportions
  prop <- .checkPropInLearner(labels, prop)
  
  # check nbIter
  if( nbIter < 1) { stop("nbIter must be greater or equal to 1")}
  
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?learnGamma for the list of valid criterion")}
  
  # check models
  if (!clusterValidGammaNames(models))
  { stop("models is not valid. See ?clusterGammaNames for the list of valid model names")}
  # Create model
  model = .createMixtureModel("ClusterGamma", data, labels, prop);
  # Create algorithm
  algo = learnAlgo( algo, nbIter, epsilon)
  # start estimation of the models
  resFlag = .Call("learnMixture", model, models, algo, criterion, nbCore, PACKAGE="MixAll");
  # set names
  if (resFlag != TRUE ) {cat("WARNING: An error occur during the learning process");}
  colnames(model@component@shape) <- colnames(model@component@data);
  colnames(model@component@scale) <- colnames(model@component@data);
  model
}

#' @rdname learners
learnCategorical <- function( data, labels, prop = NULL
                            , models=clusterCategoricalNames(prop = "equal")
                            , algo="simul", nbIter = 100, epsilon = 1e-08
                            , criterion="ICL"
                            , nbCore = 1)
{
  # check data
  labels = .checkDataInLearner(data, labels)
  # check proportions
  prop <- .checkPropInLearner(labels, prop)
  
  # check nbIter
  if( nbIter < 1) { stop("nbIter must be greater or equal to 1")}
  
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?learnGamma for the list of valid criterion")}
  
  # check models
  if (!clusterValidCategoricalNames(models))
  { stop("models is not valid. See ?clusterCategoricalNames for the list of valid model names")}
  # Create model
  model = .createMixtureModel("ClusterCategorical", data, labels, prop);
  # Create algorithm
  algo = learnAlgo( algo, nbIter, epsilon)
  # start estimation of the models
  resFlag = .Call("learnMixture", model, models, algo, criterion, nbCore, PACKAGE="MixAll");
  # set names
  if (resFlag != TRUE ) {cat("WARNING: An error occur during the learning process");}
  dim(model@component@plkj) <- c(model@component@nbModalities, model@nbCluster, ncol(data));
  model
}

#' @param dim integer giving the dimension of the Gaussian density. Default is 10.
#' @param kernelName string with a kernel name. Possible values:
#' "gaussian", "polynomial", "exponential", "linear", "hamming". Default is "gaussian".
#' @param kernelParameters [\code{\link{vector}}] with the parameters of
#' the chosen kernel. Default is 1.
#' 
#' @rdname learners

learnKernel <- function( data, labels, prop = NULL
                       , models=clusterKernelNames(prop = "equal")
                       , algo="impute", nbIter = 100, epsilon = 1e-08
                       , dim = 10, kernelName = "gaussian", kernelParameters = 1.
                       , criterion="ICL"
                       , nbCore = 1)
{
  # check data
  labels = .checkDataInLearner(data, labels)
  # check proportions
  prop <- .checkPropInLearner(labels, prop)
  
  # check nbIter
  if( nbIter < 1) { stop("nbIter must be greater or equal to 1")}
  
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?learnGamma for the list of valid criterion")}
 
  # check models
  if (!clusterValidKernelNames(models))
  { stop("models is not valid. See ?clusterKernelNames for the list of valid model names")}
  
  # check kernelName
  if(sum(kernelName %in% c("gaussian","polynomial", "exponential","linear","hamming")) != 1)
  { stop("kernelName is not valid. See ?clusterKernelNames for the list of valid kernel name")}
  if (is.null(kernelParameters)) { kernelParameters = c(1)}
  
  # Create model
  model = new("ClusterKernel", data, dim, kernelName, kernelParameters)
  # update fields
  tik <- matrix(0, nrow=length(labels), ncol=length(levels(labels)))
  for (i in 1:length(labels))
  { tik[i, as.numeric(labels[i])] <- 1}
  model@zi       = as.integer(labels) -1L; # base 0 for the labels class
  model@tik      = tik;
  model@pk       = prop
  # Create algorithm
  algo = learnAlgo( algo, nbIter, epsilon)
  # start estimation of the models
  resFlag = .Call("learnKernelMixture", model, models, algo, criterion, nbCore, PACKAGE="MixAll");
  # set names
  if (resFlag != TRUE ) {cat("WARNING: An error occur during the learning process");}
  model
}

#-----------------------------------------------------------------------
#' This function learn the optimal mixture model when the class labels are known
#' according to the \code{criterion} among the list of model given in \code{models}.
#'
#' @param data [\code{list}] containing the data sets (matrices and/or data.frames).
#' If data sets contain NA values, these missing values will be estimated during
#' the estimation process.
#' @param models either a [\code{vector}] of character or a [\code{list}] of
#' same length than data. If \code{models} is a vector, it contains the model
#' names to use in order to fit each data set. If \code{models} is a list, it
#' must be of the form 
#' \code{models = list( modelName, dim, kernelName, modelParameters) }
#' Only modelName is required.
#' @param labels vector or factors giving the label class.
#' @param prop [\code{\link{vector}}] with the proportions of each class.
#' If NULL the proportions will be estimated using the labels.
#' @param algo character defining the algo to used in order to learn the model.
#' Possible values: "simul" (default), "impute" (faster but can produce biased results).
#' @param nbIter integer giving the number of iterations to do.
#' algo is "impute" this is the maximal authorized number of iterations. Default is 100.
#' @param epsilon real giving the variation of the log-likelihood for stopping the
#' iterations. Not used if algo is "simul". Default value is 1e-08.
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL" (default).
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the heart disease data set
#' data(HeartDisease.cat)
#' data(HeartDisease.cont)
#' ## with default values
#' ldata = list(HeartDisease.cat, HeartDisease.cont);
#' models = c("categorical_pk_pjk","gaussian_pk_sjk")
#' model <- clusterMixedData(ldata, models, nbCluster=2:5, strategy = clusterFastStrategy())
#'
#' ## get summary
#' summary(model)
#'
#' ## get estimated missing values
#' missingValues(model)
#'
#' \dontrun{
#' ## print model
#' print(model)
#' ## use graphics functions
#' plot(model)
#' }
#'
#' @return An instance of the [\code{\linkS4class{ClusterMixedData}}] class.
#' @author Serge Iovleff
#' @export
#'
learnMixedData <- function( data, models, labels, prop = NULL
                          , algo="impute", nbIter = 100, epsilon = 1e-08
                          , criterion="ICL"
                          , nbCore = 1)
{
  # check data and models
  if (!is.list(data))     { stop("data must be a list");}
  if (!is.vector(models)) { stop("models must be a vector of character");}
  if (length(data) != length(models)) { stop("data and models must be of equal lengths");}
  
  labels = as.factor(as.vector(labels))
  prop<- .checkPropInLearner(labels, prop)
  nbCluster <-length(prop)
  # create list of component
  ldata <- vector("list", length(data));
  for (i in 1:length(data))
  {
    # get the name of the model and the parameters for the kernel if it is a list
    if (is.list(models))
    {
      param <- models[[i]];
      if (is.list(param)) { modelName <- param$modelName}
      else {modelName <- param}
    }
    else # otherwise param and modelName are the same
    {
      param <- models[i];
      modelName <- models[i];
    }
    # check if it is a Categorical model 
    if( clusterValidCategoricalNames(modelName) )
    { ldata[[i]] <- new("ClusterCategoricalComponent", data[[i]], nbCluster, modelName);}
    else 
    {  # check if it is a Gamma model
      if( clusterValidGammaNames(modelName) )
      { ldata[[i]] <- new("ClusterGammaComponent", data[[i]], nbCluster, modelName);}
      else
      { # check if it is a diagonal Gaussian model
        if( clusterValidDiagGaussianNames(modelName) )
        { ldata[[i]] <- new("ClusterDiagGaussianComponent", data[[i]], nbCluster, modelName);}
        else
        { # check if it a kernel model
          if( clusterValidKernelNames(modelName) )
          {
            if (is.list(param)) # get kernel parameters and check values
            {
              kernelName = param$kernelName;
              if(is.null(kernelName)) { kernelName = "gaussian";}
              kernelParameters = param$kernelParameters;
              if(is.null(kernelParameters)) { kernelParameters = c();}
              dim = param$dim;
              if(is.null(dim)) { dim = 10;}
            }
            else # set default values
            {
              kernelName <- "gaussian"
              kernelParameters <- c()
              dim <- 10
            }
            # check kernel name
            if(sum(kernelName %in% c("gaussian","polynomial", "exponential","linear","hamming")) != 1)
            { stop("kernelName is not valid. See ?clusterKernel for the list of valid kernel name")}
            # check dim
            if (dim < 1) { stop("The dimension must be greater or equal to 1")}
            # create component
            component <- new("ClusterKernelComponent", data[[i]], dim= dim, nbCluster= nbCluster, modelName= modelName);
            # compute gram matrix
            resFlag  <- .Call("clusterKernelCompute", component, kernelName, kernelParameters, PACKAGE="MixAll");
            if (!resFlag) { stop("error in gram matrix computation");}
            # restore original data set and update ldata list
            component@rawData <- data[[i]];
            ldata[[i]] <- component;
          }
          else
          {
            stop("invalid model name");
          }
        } # else diag Gaussian
      } # else gamma
    } # else categorical
  } # for i
  # create model
  model = .createMixtureModel("ClusterMixedData", ldata, labels, prop);
  # Create algorithm
  algo = learnAlgo( algo, nbIter, epsilon)
  # start estimation of the models
  resFlag  <- FALSE;
  if (length(nbCluster) >0)
  {
    resFlag = .Call("learnMixedData", model, algo, criterion, nbCore, PACKAGE="MixAll");
  }
  # set names
  if (resFlag != TRUE) {cat("WARNING: An error occurs during the clustering process");}
  for (i in 1:length(data))
  {
    if(clusterValidCategoricalNames(models[i]))
    { dim(model@ldata[[i]]@plkj) <- c(model@ldata[[i]]@nbModalities, model@nbCluster, ncol(model@ldata[[i]]@data))}
  }
  model
}

.createMixtureModel <- function(name, data, labels, prop)
{
  # create model
  model = new(Class = name, data, nbCluster = nlevels(labels))
  # update fields
  tik <- matrix(0, nrow=length(labels), ncol=length(levels(labels)))
  for (i in 1:length(labels))
  { tik[i, as.numeric(labels[i])] <- 1}
  model@zi       = as.integer(labels) -1L; # base 0 for the labels class
  model@tik      = tik;
  model@pk       = prop
  model
}

.checkDataInLearner <- function(data, labels)
{
  data   = as.matrix(data)
  labels = as.vector(labels)
  if (ncol(data) < 1) {stop("Error: empty data set")}
  if (nrow(data) != length(labels)) {stop("Error: data and labels must have the same lengths")}
  as.factor(labels)
}

.checkPropInLearner <- function(labels, prop)
{
  # check labels and proportions
  if (is.null(prop))
  { prop<-as.vector(table(labels))}
  else
  {
    prop <- as.vector(prop);
    if (nlevels(labels) != length(prop))
    { stop("the number of levels in labels does not match the proportions length")}
  }
  if (min(prop) <= 0) stop("All proportions must be strictly greater than 0");
  prop/sum(prop);
}
