#-----------------------------------------------------------------------
#     Copyright (C) 2012-2016  Inria
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
#' @include ClusterStrategy.R
NULL

#' Interface base Class [\code{\linkS4class{IClusterModel}}] for Cluster models.
#'
#' This class encapsulate the common parameters of all the Cluster models.
#'
#' A Cluster model is a model of the form
#' \deqn{
#'   f({x}|\boldsymbol{\theta})
#'     \sum_{k=1}^K p_k h({x};\boldsymbol{\lambda}_k,\boldsymbol{\alpha})
#'    \quad {x} \in J.
#' }
#' where h can be either a pdf, a discrete probability, (homogeneous case)
#' or a product of arbitrary pdf and discrete probabilities (mixed data case).
#'
#' @slot nbSample  Integer with the number of samples of the model.
#' @slot nbCluster Integer with the number of cluster of the model.
#' @slot pk        Vector of size K with the proportions of each mixture.
#' @slot tik       Matrix of size \eqn{n \times K} with the posterior probability of
#' the ith individual to belong to kth cluster.
#' @slot lnFi        Vector of size n with the log-likelihood of the ith individuals.
#' @slot zi        Vector of integer of size n  with the attributed class label of the individuals.
#' @slot lnLikelihood Real given the ln-liklihood of the Cluster model.
#' @slot criterion Real given the value of the AIC, BIC or ICL criterion.
#' @slot nbFreeParameter Integer given the number of free parameters of the model.
#' @slot strategy  the instance of the [\code{\linkS4class{ClusterStrategy}}] used in the
#' estimation process of the mixture. Default is clusterStrategy().
#'
#' @examples
#'   getSlots("IClusterModel")
#'
#' @author Serge Iovleff
#'
#' @name IClusterModel
#' @rdname ClusterModels-class
#' @aliases IClusterModel-class
#' @exportClass IClusterModel
setClass(
  Class="IClusterModel",
  representation( nbSample = "numeric"
                , nbCluster = "numeric"
                , pk = "numeric"
                , tik = "matrix"
                , lnFi = "numeric"
                , zi = "integer"
                , lnLikelihood = "numeric"
                , criterion = "numeric"
                , nbFreeParameter = "numeric"
                , strategy = "ClusterStrategy"
                , "VIRTUAL"
                ),
  # validity function
  validity=function(object)
  {
    nbSample  = object@nbSample
    nbCluster = object@nbCluster

    # check nbSample
    if (round(nbSample)!=object@nbSample)
    {stop("nbSample must be an integer.")}
    if( nbSample < 2 )
    { stop("nbSample must be greater than 1.")}

    # check nbCluster
    if (round(nbCluster)!=object@nbCluster)
    {stop("nbCluster must be an integer.")}
    if( nbCluster < 1 )
    { stop("nbCluster must be greater than 0.")}

    # check pk
    if (length(object@pk) != nbCluster)
    {stop("pk must have length nbCluster.")}

    # check tik
    if (ncol(object@tik) != nbCluster)
    {stop("tik must have nbCluster columns.")}
    if (nrow(object@tik) != nbSample)
    {stop("tik must have nbSample rows.")}

    # check lnFi
    if (length(object@lnFi) != nbSample)
    {stop("fi must have nbSample size.")}

    # check zi
    if (length(object@zi) != nbSample)
    {stop("zi must have nbSample size.")}
    return(TRUE)

    # check nbFreeParameter
    if (round(object@nbFreeParameter)!=object@nbFreeParameter)
    {stop("nbFreeParameter must be an integer.")}

    if(class(strategy)[1] != "ClusterStrategy")
    {stop("strategy is not a Cluster Strategy class (must be an instance of the class ClusterStrategy).")}

    return(TRUE)

  }
)

#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{IClusterModel}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
#'
setMethod(
  f="initialize",
  signature=c("IClusterModel"),
  definition=function(.Object, nbSample, nbCluster)
  {
    # for nbSample
    if(missing(nbSample)) { stop("nbSample is mandatory in IClusterModel.")}
    .Object@nbSample <- nbSample;
    # for nbCluster
    if(missing(nbCluster)) { stop("nbCluster is mandatory in IClusterModel.")}
    .Object@nbCluster<-nbCluster
    # resize
    .Object@pk   <- rep(1/nbCluster, nbCluster)
    .Object@tik  <- matrix(1/nbCluster, nbSample, nbCluster)
    .Object@lnFi <- rep(0, nbSample)
    .Object@zi   <- as.integer(rep(1, nbSample))
    # set default values
    .Object@lnLikelihood = -Inf
    .Object@criterion    =  Inf
    .Object@nbFreeParameter = 0
    # set strategy
    .Object@strategy <- clusterStrategy()
    # validObject(.Object) will be called at the end of the initialization process
    # in the derived classes
    return(.Object)
  }
)

#' @rdname print-methods
#' @aliases print print,IClusterModel-method
#'
setMethod(
  f="print",
  signature=c("IClusterModel"),
  function(x,...)
  {
    cat("* nbSample       = ", x@nbSample, "\n")
    cat("* nbCluster      = ", x@nbCluster, "\n")
    cat("* lnLikelihood   = ", x@lnLikelihood,"\n")
    cat("* nbFreeParameter= ", x@nbFreeParameter,"\n")
    cat("* criterion      = ", x@criterion, "\n")
    cat("* zi             =\n")
    print( format(x@zi), quote=FALSE)
  }
)

#' @rdname show-methods
#' @aliases show show,IClusterModel-method
setMethod(
  f="show",
  signature=c("IClusterModel"),
  function(object)
  {
    cat("* nbSample       = ", object@nbSample, "\n")
    cat("* nbCluster      = ", object@nbCluster, "\n")
    cat("* lnLikelihood   = ", object@lnLikelihood,"\n")
    cat("* nbFreeParameter= ", object@nbFreeParameter,"\n")
    cat("* criterion      = ", object@criterion, "\n")
  }
)

#' @rdname summary-methods
#' @aliases summary summary,IClusterModel-method
setMethod(
  f="summary",
  signature=c("IClusterModel"),
  function(object,...)
  {
    cat("* nbSample       = ", object@nbSample, "\n")
    cat("* nbCluster      = ", object@nbCluster, "\n")
    cat("* lnLikelihood   = ", object@lnLikelihood,"\n")
    cat("* nbFreeParameter= ", object@nbFreeParameter, "\n")
    cat("* criterion      = ", object@criterion, "\n")
  }
)

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{IClusterComponent}}] class
#'
#' This class defines a component of a mixture Model
#'
#' @slot data  Matrix with the data set
#' @slot missing  Matrix with the indexes of the missing values
#' @slot modelName model name associated with the data set
#'
#' @examples
#' getSlots("IClusterComponent")
#'
#' @author Serge Iovleff
#'
#' @name IClusterComponent
#' @rdname IClusterComponent-class
#' @aliases IClusterComponent-class
#' @exportClass IClusterComponent
#'
setClass(
  Class="IClusterComponent",
  representation( data = "matrix"
                , missing = "matrix"
                , modelName = "character"
                , "VIRTUAL"),
  validity=function(object)
  {
# called too soon when Component is part of an other S4 class
#    if (!is.matrix(data)) { stop("data must be a matrix in IClusterComponent");}
#    if (!is.matrix(missing)) { stop("missing must be a matrix in IClusterComponent");}
#    if (!is.character(modelName)) { stop("modelName must be a character in IClusterComponent");}
    return(TRUE);
  }
)
#' Initialize an instance of a MixAll S4 class.
#'
#' Initialization method of the [\code{\linkS4class{IClusterComponent}}] class.
#' Used internally in the 'MixAll' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
  f="initialize",
  signature=c("IClusterComponent"),
  definition=function(.Object, data, modelName)
  {
    # fill data missing and modelName
    .Object@data      <- as.matrix(data);
    .Object@missing   <- which(is.na(.Object@data), arr.ind=TRUE);
    .Object@modelName <- modelName;
    return(.Object)
  }
)
#' @rdname print-methods
#' @aliases print print,IClusterComponent-method
#'
setMethod(
  f="print",
  signature=c("IClusterComponent"),
  function(x,...)
  {
    cat("* model name   =",x@modelName,"\n");
    cat("* data         =\n");
    print(format(x@data),quote=FALSE)
    cat("* missing     =\n");
    print(format(x@missing),quote=FALSE)
  }
)

#' @rdname show-methods
#' @aliases show-IClusterComponent,IClusterComponent,IClusterComponent-method
setMethod(
  f="show",
  signature=c("IClusterComponent"),
  function(object)
  {
    cat("* model name   =",object@modelName,"\n");
    cat("*\n");
    if(length(object@data)!=0)
    {
      nrowShow <- min(10,nrow(object@data));
      ncolShow <- min(10,ncol(object@data));
      cat("* data (limited to 10 samples and 10 variables) =\n")
      print(format(object@data[1:nrowShow,1:ncolShow]),quote=FALSE)
    }
    cat("* ... ...\n")
  }
)

#' @rdname summary-methods
#' @aliases summary summary,IClusterComponent-method
setMethod(
  f="summary",
  signature=c("IClusterComponent"),
  function(object, ...)
  { cat("* model name     =",object@modelName,"\n");}
)
## 
## #-----------------------------------------------------------------------
## #' Interface base Class [\code{\linkS4class{IClusterPredict}}] for predictors
## #' using cluster models.
## #'
## #' This class encapsulate the common parameters of all the predictors using a
## #' cluster model.
## #'
## #'
## #' @slot tik  Matrix of size \eqn{n \times K} with the predicted posterior
## #' probabilities of the ith individual to belong to kth cluster.
## #' @slot lnFi Vector of size n with the predicted log-likelihood of the ith
## #' individuals.
## #' @slot zi   Vector of integer of size n with the predicted class label of
## #' the individuals.
## #'
## #' @examples
## #'   getSlots("IClusterPredict")
## #'
## #' @author Serge Iovleff
## #'
## #' @name IClusterPredict
## #' @rdname IClusterPredict-class
## #' @aliases IClusterPredict-class
## #' @exportClass IClusterPredict
## #'
## setClass(
##     Class="IClusterPredict",
##     representation( nbSample = "numeric"
##                   , nbCluster = "numeric"
##                   , tik = "matrix"
##                   , lnFi = "numeric"
##                   , zi = "integer"
##                   , "VIRTUAL"
##                   ),
##     validity=function(object)
##     {
##       return(TRUE)
##     }
## )
## 
## #' Initialize an instance of a MixAll S4 class.
## #'
## #' Initialization method of the [\code{\linkS4class{IClusterPredict}}] class.
## #' Used internally in the 'MixAll' package.
## #'
## #' @rdname initialize-methods
## #' @keywords internal
## #'
## setMethod(
##     f="initialize",
##     signature=c("IClusterPredict"),
##     definition=function(.Object, nbSample, nbCluster)
##     {
##       # for nbSample
##       if(missing(nbSample)) { stop("nbSample is mandatory in IClusterPredict.")}
##       .Object@nbSample <- nbSample;
##       # for nbCluster
##       if(missing(nbCluster)) { stop("nbCluster is mandatory in IClusterPredict.")}
##       .Object@nbCluster<-nbCluster
##       # resize
##       .Object@pk   <- rep(1/nbCluster, nbCluster)
##       .Object@tik  <- matrix(1/nbCluster, nbSample, nbCluster)
##       .Object@lnFi <- rep(0, nbSample)
##       .Object@zi   <- as.integer(rep(1, nbSample))
##       # validObject(.Object) will be called at the end of the initialization process
##       # in the derived classes
##       return(.Object)
##     }
## )
## 
