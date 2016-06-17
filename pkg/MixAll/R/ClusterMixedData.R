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
#' @include ClusterModelNames.R IClusterModel.R
NULL

#-----------------------------------------------------------------------
#' Create an instance of the [\code{\linkS4class{ClusterMixedData}}] class
#'
#' This function computes the optimal mixture model for mixed data according
#' to the \code{criterion} among the number of clusters given in
#' \code{nbCluster} using the strategy specified in [\code{strategy}].
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
#' @param nbCluster [\code{\link{vector}}] with the number of clusters to test.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing
#' the strategy to run. Default is clusterStrategy().
#' @param criterion character defining the criterion to select the best model.
#' The best model is the one with the lowest criterion value.
#' Possible values: "BIC", "AIC", "ICL". Default is "ICL".
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
clusterMixedData <- function( data, models, nbCluster=2
                            , strategy=clusterFastStrategy()
                            , criterion="ICL"
                            , nbCore = 1)
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin   = min(nbCluster);
  nbClusterMax   = max(nbCluster);
  if (nbClusterMin < 1) { stop("The number of clusters must be greater or equal to 1")}
  # check criterion
  if(sum(criterion %in% c("BIC","AIC", "ICL")) != 1)
  { stop("criterion is not valid. See ?clusterMixedData for the list of valid criterion")}
  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);
  # check data and models
  if (!is.list(data))     { stop("data must be a list");}
  if (!is.vector(models)) { stop("models must be a vector of character");}
  if (length(data) != length(models)) { stop("data and models must be of equal lengths");}

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
    { ldata[[i]] <- new("ClusterCategoricalComponent", data[[i]], nbClusterMin, modelName);}
    else 
    {  # check if it is a Gamma model
      if( clusterValidGammaNames(modelName) )
      { ldata[[i]] <- new("ClusterGammaComponent", data[[i]], nbClusterMin, modelName);}
      else
      { # check if it is a diagonal Gaussian model
        if( clusterValidDiagGaussianNames(modelName) )
        { ldata[[i]] <- new("ClusterDiagGaussianComponent", data[[i]], nbClusterMin, modelName);}
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
            component <- new("ClusterKernelComponent", data[[i]], dim= dim, nbCluster= nbClusterMin, modelName= modelName);
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
  # Create model
  model = new("ClusterMixedData", ldata)
  model@strategy = strategy;
  # start estimation of the models
  resFlag  <- FALSE;
  if (length(nbCluster) >0)
  {
   resFlag = .Call("clusterMixedData", model, nbCluster, strategy, criterion, nbCore, PACKAGE="MixAll");
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

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterMixedData}}] class
#'
#' This class defines a mixed data mixture Model.
#'
#' This class inherits from the [\code{\linkS4class{IClusterModel}}] class.
#' A model for mixed data is a mixture model of the form:
#' \deqn{
#' f({{x}}_i=({{x}}_{1i}, {{x}}_{2i},\ldots {{x}}_{Li})|\theta)
#' = \sum_{k=1}^K p_k \prod_{l=1}^L h({{x}}_{li}| \lambda_{lk},\alpha_l).
#' }
#' The density functions (or probability distribution functions)
#' \deqn{h(.|\lambda_{lk},\alpha_l)}
#' can be any implemented model (Gaussian, Poisson,...).
#'
#' @slot ldata  a list of IClusterComponent.
#' @seealso [\code{\linkS4class{IClusterModel}}] class
#'
#' @examples
#' getSlots("ClusterMixedData")
#'
#' @author Serge Iovleff
#'
#' @name ClusterMixedData
#' @rdname ClusterMixedData-class
#' @aliases ClusterMixedData-class
#' @exportClass ClusterMixedData
#'
setClass(
    Class="ClusterMixedData",
    representation( ldata = "list"),
    contains=c("IClusterModel"),
    validity=function(object)
    {
      nbData = length(object@ldata)
      if (nbData == 0) {stop("At least on data set must be given.");}
      for (l in 1:nbData)
      {
        if (nrow(object@ldata[[1]]@data) != object@nbSample)
        {stop("All data sets must have the same number of individuals (number of rows).");}
      }
      return(TRUE)
    }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterMixedData}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterMixedData"),
    definition=function(.Object, ldata =list(), nbCluster=2)
    {
      # for data
      if(missing(ldata)) {stop("ldata is mandatory in ClusterMixedData.")}
      nbData = length(ldata)
      if (nbData == 0) {stop("At least on data set must be given.")}
      .Object@ldata <- ldata;
      # take first element of the list, this will give us the dimensions
      nbSample = nrow(.Object@ldata[[1]]@data);
      .Object <- callNextMethod(.Object, nbSample, nbCluster)
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterMixedData-method
#'
setMethod(
  f="print",
  signature=c("ClusterMixedData"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod()
    nbData <- length(x@ldata)
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        cat("* model name = ", x@ldata[[l]]@modelName, "\n")
        print(format(x@ldata[[l]]@data),quote=FALSE);
      }
    }
    cat("****************************************\n")
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        for(k in 1:length(x@pk))
        {
          cat("*** Cluster: ",k,"\n")
          cat("* Proportion = ", format(x@pk[k]), "\n")
          cat("* model name = ", x@ldata[[l]]@modelName, "\n");
          print(x@ldata[[l]],k);
        }
      }
      cat("****************************************\n")
    }
  }
)

#' @rdname show-methods
#' @aliases show-ClusterMixedData,ClusterMixedData,ClusterMixedData-method
setMethod(
    f="show",
    signature=c("ClusterMixedData"),
    function(object)
    {
      cat("****************************************\n")
      callNextMethod()
      nbData <- length(object@ldata)
      if(nbData>0)
      {
        for (l in 1:nbData)
        {
          cat("* model name = ", object@ldata[[l]]@modelName, "\n")
          nrowShow <- min(10,nrow(object@ldata[[l]]@data))
          ncolShow <- min(10,ncol(object@ldata[[l]]@data))
          cat("* data (limited to 10 samples and 10 variables) =\n")
          print(format(object@ldata[[l]]@data[1:nrowShow,1:ncolShow]),quote=FALSE)
        }
      }
      cat("* ... ...\n")

      cat("****************************************\n")
      if(nbData>0)
      {
        for (l in 1:nbData)
        {
          for(k in 1:length(object@pk))
          {
            cat("*** Cluster: ",k,"\n")
            cat("* Proportion = ", format(object@pk[k]),"\n")
            cat("*\n")
            cat("* model name = ", object@ldata[[l]]@modelName, "\n");
            print(object@ldata[[l]], k);
          }
        }
        cat("****************************************\n")
      }
    }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterMixedData-method
setMethod(
    f="summary",
    signature=c("ClusterMixedData"),
    function(object, ...)
    {
      cat("**************************************************************\n")
      callNextMethod()
      cat("**************************************************************\n")
    }
)

#' Plotting of a class [\code{\linkS4class{ClusterMixedData}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterMixedData}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterMixedData}}]
#' @param y a number between 1 and K-1.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @aliases plot-ClusterMixedData
#' @docType methods
#' @rdname plot-ClusterMixedData-method
#' @export
#'
#' @seealso \code{\link{plot}}
#' @examples
#' \dontrun{
#'   ## the car data set
#'   data(car)
#'   model <- clusterMixedData(car, 3, strategy = clusterFastStrategy())
#'   plot(model)
#'   }
#'
setMethod(
    f="plot",
    signature=c("ClusterMixedData"),
    function(x, y, ...)
    {
      # total number of cluster in the data set
      nbCluster = ncol(x@tik);
      # check y, no y => display all dimensions
      if (missing(y)) { y=1:(nbCluster-1); }
      else
      { if (round(y)!=y) {stop("y must be an integer.")}
        if (y>nbCluster-1)
          stop("y should not be greater than K-1")
        y <- 1:y
      }
      # get representation
      Y=.visut(x@tik, nbCluster);
      if (nbCluster == 2) { ndim = 1;}
      else { ndim = ncol(Y);}
      # Compute gaussian statistics
      mean  <- matrix(0, nrow = x@nbCluster, ncol =ndim)
      sigma <- matrix(1, nrow = x@nbCluster, ncol =ndim)
      for (k in 1:nbCluster)
      {
        wcov = cov.wt(as.matrix(Y), x@tik[,k], method = "ML");
        mean[k,]  = wcov$center;
        sigma[k,] = sqrt(diag(wcov$cov))
      }
      # create gaussian model
      gauss<-new("ClusterDiagGaussian", Y, nbCluster = x@nbCluster)
      gauss@component@mean = mean
      gauss@component@sigma= sigma
      gauss@pk   = x@pk
      gauss@tik  = x@tik
      gauss@lnFi = x@lnFi
      gauss@zi   = x@zi
      #gauss@component@missing     = x@component@missing
      gauss@lnLikelihood = x@lnLikelihood
      gauss@criterion    = x@criterion
      gauss@nbFreeParameter = x@nbFreeParameter
      gauss@strategy        = x@strategy
      .clusterPlot(gauss, y, .dGauss,...);
    }
)

# get logistic representation
.visut <- function(t, gp)
{ m <- min(t[,gp]);
  if (m==0) t[,gp] = t[,gp] + 1e-30;
  return(scale(log(sweep(t,1,t[,gp],FUN="/")+ 1e-30), center=TRUE, scale=FALSE)[,-gp])
}

