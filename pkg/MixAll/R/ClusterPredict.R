#-----------------------------------------------------------------------
#     Copyright (C) 2012-2018  Serge Iovleff, University Lille 1, Inria
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
NULL

#-----------------------------------------------------------------------
#' Create an instance of a mixture model S4 class with predicted values
#'
#' This function predicts the best cluster each sample in data belongs to. 
#'
#' @param data frame or matrix containing the data. Rows correspond to observations
#' and columns correspond to variables. If the data set contains NA values, they
#' will be estimated during the predicting process.
#' @param model estimated clustering model
#' @param algo an instance of \code{\link{ClusterAlgoPredict}} S4 class
#' @param nbCore integer defining the number of processors to use (default is 1, 0 for all).
#'
#' @examples
#' ## A quantitative example with the famous iris data set
#' data(iris)
#' ## add 10 missing values as random
#' x = as.matrix(iris[1:4]); n <- nrow(x); p <- ncol(x);
#' indexes <- matrix(c(round(runif(5,1,n)), round(runif(5,1,p))), ncol=2);
#' x[indexes] <- NA;
#' ## sample train and test data sets
#' indexes <- sample(1:nrow(x), nrow(x)/2)
#' train <- x[ indexes,]
#' test  <- x[-indexes,]
#' ## estimate model (using fast strategy, results may be misleading)
#' model1 <- clusterDiagGaussian( data =train, nbCluster=2:3
#'                              , models=c( "gaussian_pk_sjk")
#'                              , strategy = clusterFastStrategy()
#'                              )
#' ## get summary
#' summary(model1)
#' ## compute prediction and compare
#' model2 <- clusterPredict(test, model1)
#' show(model2)
#' as.integer(iris$Species[-indexes])
#'
#' @return An instance of [\code{\linkS4class{ClusterPredict}}] with predicted
#' values
#' @author Serge Iovleff
#'
#'
clusterPredict <- function( data, model, algo = clusterAlgoPredict(), nbCore = 1)
{
  # for data
  if(missing(data)) { stop("data is mandatory in clusterPredict.")}
  data = as.matrix(data)
  if (nrow(data)<1) { stop("data is empty in clusterPredict.")}
  
  # for model
  if(missing(model)) { stop("model is mandatory in clusterPredict.")}
  if(!is(model,"IClusterModel")) { stop("model must be an instance of IClusterModel or a derived class.")}
  
  # check dimension
  if (ncol(model@component@data) != ncol(data))
  { stop("data does not have the same number of variables than the model;")}
  
  result = new("ClusterPredict", data, model@nbCluster, algo)
  # nothing to do
  if (nrow(data)>=1)
  {
    # start estimation of the models
    resFlag = .Call("clusterPredict", model, result, PACKAGE="MixAll");
    # set names
    if (resFlag != TRUE ) {cat("WARNING: An error occur during the clustering process");}
  }
  result
}


#' Class [\code{\linkS4class{ClusterPredict}}] for predicting 
#'
#' This class encapsulate the parameters for predicted data.
#'
#' @slot data  Matrix with the data set
#' @slot missing   Matrix with the indexes of the missing values
#' @slot nbSample  Integer with the number of samples 
#' @slot nbCluster Integer with the number of cluster
#' @slot pk        Vector of size K with the proportions of each mixture.
#' @slot tik       Matrix of size \eqn{n \times K} with the posterior probability of
#' the ith individual to belong to kth cluster.
#' @slot lnFi      Vector of size n with the log-likelihood of the ith individuals.
#' @slot zi        Vector of integer of size n  with the attributed class label of the individuals
#' @slot algo      an instance of [\code{\linkS4class{ClusterAlgoPredict}}] 
#'
#'
#' @examples
#'   getSlots("ClusterPredict")
#'
#' @author Serge Iovleff
#'
#' @name ClusterPredict
#' @rdname ClusterPredict-class
#' @aliases ClusterPredict-class
#' @exportClass ClusterPredict
setClass(
  Class="ClusterPredict",
  # members
  representation( data = "matrix"
                , missing   = "matrix"
                , nbSample  = "numeric"
                , nbCluster = "numeric"
                , pk        = "numeric"
                , tik       = "matrix"
                , lnFi      = "numeric"
                , zi        = "integer"
                , algo      = "ClusterAlgoPredict"
              ),
  # validity function
  validity=function(object)
  {
    nbSample  <- object@nbSample
    nbCluster <- object@nbCluster
    
    # check nbSample
    if (round(nbSample) != nbSample)
    {stop("Error in ClusterPredict validity. nbSample must be an integer.")}
    
    # check nbCluster
    if (round(nbCluster)!= nbCluster)
    {stop("Error in ClusterPredict validity. nbCluster must be an integer.")}
    if( nbCluster < 1 )
    { stop("Error in ClusterPredict validity. nbCluster must be greater than 0.")}
    
    # check data
    if (nrow(object@data) != nbSample)
    {stop("Error in ClusterPredict validity. data must have nbSample rows.")}
    
    # check pk
    if (length(object@pk) != nbCluster)
    { stop("Error in ClusterPredict validity. pk must have length nbCluster.")}
    
    # check tik
    if (ncol(object@tik) != nbCluster)
    { stop("Error in ClusterPredict validity. tik must have nbCluster columns.")}
    if (nrow(object@tik) != nbSample)
    { stop("Error in ClusterPredict validity. tik must have nbSample rows.")}
    
    # check lnFi
    if (length(object@lnFi) != nbSample)
    { stop("Error in ClusterPredict validity. lnFi must have nbSample size.")}
    
    # check zi
    if (length(object@zi) != nbSample)
    { stop("Error in ClusterPredict validity. zi must have nbSample size.")}
    
    # check algo
    if (class(object@algo)[1] != "ClusterAlgoPredict")
    { stop("Error in ClusterPredict validity. algo must be an instance of ClusterPredictAlgo.")}
    return(TRUE)
  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterPredict}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
#'
setMethod(
  f="initialize",
  signature=c("ClusterPredict"),
  definition=function(.Object, data, nbCluster, algo = clusterAlgoPredict())
  {
    # for data
    if(missing(data)) { stop("data is mandatory in ClusterPredict.")}
    .Object@data <- data
    .Object@nbSample <- nrow(data)
    .Object@missing  <- which(is.na(.Object@data), arr.ind=TRUE);

    # for nbCluster
    if(missing(nbCluster)) { stop("nbCluster is mandatory in ClusterPredict.")}
    .Object@nbCluster<-nbCluster
    
    # create arrays
    .Object@pk   <- rep(1/nbCluster, nbCluster)
    .Object@tik  <- matrix(1/nbCluster, .Object@nbSample, nbCluster)
    .Object@lnFi <- rep(0, .Object@nbSample)
    .Object@zi   <- as.integer(rep(1, .Object@nbSample))
    .Object@algo <- algo
    
    # valid object
    validObject(.Object)
    
    # in the derived classes
    return(.Object)
  }
)

#' @rdname print-methods
#' @aliases print print,ClusterPredict-method
#'
setMethod(
  f="print",
  signature=c("ClusterPredict"),
  function(x,...)
  {
    cat("****************************************\n")
    cat("* nbSample       = ", x@nbSample, "\n")
    cat("* nbCluster      = ", x@nbCluster, "\n")
    cat("* zi             =\n")
    print( format(x@zi), quote=FALSE)
  }
)

#' @rdname show-methods
#' @aliases show show,ClusterPredict-method
setMethod(
  f="show",
  signature=c("ClusterPredict"),
  function(object)
  {
    cat("****************************************\n")
    cat("* nbSample       = ", object@nbSample, "\n")
    cat("* nbCluster      = ", object@nbCluster, "\n")
    cat("* zi             =\n")
    print( format(object@zi), quote=FALSE)
    cat("****************************************\n")
  }
)

#' @rdname summary-methods
#' @aliases summary summary,ClusterPredict-method
setMethod(
  f="summary",
  signature=c("ClusterPredict"),
  function(object,...)
  {
    cat("****************************************\n")
    cat("* nbSample       = ", object@nbSample, "\n")
    cat("* nbCluster      = ", object@nbCluster, "\n")
    cat("****************************************\n")
  }
)


