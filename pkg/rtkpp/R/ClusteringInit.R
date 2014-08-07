#-----------------------------------------------------------------------
#     Copyright (C) 2004-2013  Serge Iovleff
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
NULL

################################################################################
#' Create an instance of [\code{\linkS4class{ClusteringInit}}] class 
#'
#' The initialization step is a two stages process: the proper initialization step
#' and some (optionnals) iterations of an algorithm [\code{\link{clusteringAlgo}}].
#'
#' @details
#' There is three ways to initialize the parameters:
#' \itemize{
#'   \item \code{random} {The initial parameters of the mixture are chosen randomly.}
#'   \item \code{class} {The initial membership of individuals are sampled randomly.}
#'   \item \code{fuzzy} {The initial probabilities of membership of individuals are
#'   sampled randomly.}
#' }
#' A few iteration of an algorithm [\code{\link{clusteringAlgo}}] are then performed.
#' It is strongly recommended to use a few number of iterations of the \code{SEM}
#' or \code{CEM} algorithms after initialization. This allow to detect "bad" initialization
#' starting point of the estimation algorithm.
#' 
#' These two stages are repeated until \code{nbInit} is reached. The initial point with the best
#' log-likelihood is conserved as the initial starting point.
#'
#' @param method Character string with the initialisation method.
#' Possible values: "random", "class", "fuzzy". Default value is "random".
#' @param algo Character string with the initialisation algorithm.
#' Possible values: "EM", "CEM", "SEM". Default value is "SEM".
#' @param nbInit integer defining the number of initialization point to test. Default value is 5.
#' @param nbIteration Integer defining the number of iteration in \code{algo}.
#' nbIteration must be a positive integer. Default values is 20. Not used if  \code{algo} = NULL.
#' @param epsilon Real defining the epsilon value for the algorithm. Default value: 0.01.
#'
#' @examples
#' clusteringInit(method = "class", nbInit=1, algo="CEM",nbIteration=50, epsilon=0.00001)
#' clusteringInit(nbIteration=0) # no algorithm
#'
#' @return a [\code{\linkS4class{ClusteringInit}}] object
#' @author Serge Iovleff
#' @export
clusteringInit <- function( method="random", nbInit=5,  algo = "SEM", nbIteration=20, epsilon=0.01)
{ return(new("ClusteringInit", nbInit=nbInit, algo=new("ClusteringAlgo", algo, nbIteration, epsilon)))}

###################################################################################
#' [\code{\linkS4class{ClusteringInit}}] class
#'
#' This class encapsulates the parameters of initialization methods of the 
#' rtkpp clustering estimation method.
#'
#' @slot method Character string with the initialization method to use. Default value: "random"
#' @slot nbInit Integer defining the number of initialization to perform. Default value: 1.
#' @slot algo An instance of \code{\linkS4class{ClusteringAlgo}} class.
#' Default value: \code{clusteringAlgo("SEM", 20, 0)}.
#'
#' @examples
#'   getSlots("ClusteringInit")
#'   new("ClusteringInit")
#'
#' @author Serge Iovleff
#' 
#' @name ClusteringInit
#' @rdname ClusteringInit-class
#' @aliases ClusteringInit-class
#' @exportClass ClusteringInit
#'
setClass(
    Class="ClusteringInit",
    slots=c(method="character", nbInit = "numeric", algo = "ClusteringAlgo"),
    prototype=list(method="random", nbInit = 5, algo = clusteringAlgo("SEM", 20, 0)),
    # validity function
    validity=function(object)
    {
      # for method
      if ( sum(object@method %in% c("random","class","fuzzy")) != 1 )
      {stop("Initialization method is not valid. See ?clusteringInit for the list of available initialization method.")}
      # for nbInit
      if (round(object@nbInit)!=object@nbInit)
      {stop("nbIInit must be an integer.")}
      if( object@nbInit < 1 ) # can't be zero
      {stop("nbInit must be strictly greater than 0.");}
      # for algo
      if (!is.null(object@algo))
      {
        if(class(object@algo)[1] != "ClusteringAlgo")
        {stop("algo is not of a clustering algorithm (must be an instance of the class ClusteringAlgo).")}
      }
      return(TRUE)
    }
)


###################################################################################
#' Create an instance of the [\code{\linkS4class{ClusteringInit}}] class using new/initialize.
#' 
#' Initialization method. Used internally in the `rtkpp' package.
#' 
# @seealso \code{\link{initialize}}
#' @keywords internal
# @rdname initialize-methods
#'
setMethod(
  f="initialize",
  signature=c("ClusteringInit"),
  definition=function(.Object,method="random",nbInit = 5,algo= clusteringAlgo("SEM", 20, 0))
  {
    # for method
    if(missing(method)) {.Object@method<-"random"}
    else  {.Object@method<-method}
    # for nbIteration
    if( missing(nbInit) ){ .Object@nbInit<-5 }
    else{.Object@nbInit<-nbInit}
    # for algo
    if(missing(algo)){ .Object@algo<-clusteringAlgo("SEM", 20, 0.1) }
    else{.Object@algo<-algo}
    # validate
    validObject(.Object)
    return(.Object)
  }
)

###################################################################################
# @name print
# @docType methods
#' @aliases print-init,ClusteringInit,ClusteringInit-method
#' @rdname print-methods
setMethod(
  f="print",
  signature=c("ClusteringInit"),
  function(x,...){
    function(object){
      cat("****************************************\n")
      cat("*** Clustering init:\n")
      cat("* method               = ", object@method, "\n")
      cat("* number of init       = ", object@nbInit, "\n")
      cat("* algorithm            = ", object@algo@algo, "\n")
      cat("* number of iterations = ", object@algo@nbIteration, "\n")
      cat("* epsilon              = ", object@algo@epsilon, "\n")
      cat("****************************************\n")
    }
  }
)

###################################################################################
# @name show
# @docType methods
#' @rdname show-methods
#' @aliases show-init,ClusteringInit,ClusteringInit-method
setMethod(
  f="show",
  signature=c("ClusteringInit"),
  function(object){
    cat("****************************************\n")
    cat("*** Clustering init:\n")
    cat("* method              = ", object@method, "\n")
    cat("* number of init      = ", object@nbInit, "\n")
    cat("* algorithm            = ", object@algo@algo, "\n")
    cat("* number of iterations = ", object@algo@nbIteration, "\n")
    cat("* epsilon              = ", object@algo@epsilon, "\n")
    cat("****************************************\n")
  }
)

###################################################################################
# @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [,ClusteringInit-method
setMethod(
  f="[", 
  signature(x = "ClusteringInit"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "method"={return(x@method)},
        "nbInit"={return(x@nbInit)},
        "algo"={return(x@algo)},
        stop("This attribute doesn't exist !")
        )
      }
    else
    {stop("This attribute is not a list !")}
  }
)

###################################################################################
#' @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [<-,ClusteringInit-method
setReplaceMethod(
  f="[", 
  signature(x = "ClusteringInit"), 
  definition=function(x,i,j,value){
    if ( missing(j) )
    {
      switch(EXPR=i,
        "method"={x@method<-value},
        "nbInit"={x@nbInit<-value},
        "algo"={x@algo<-value},
        stop("This attribute doesn't exist !")
      )
    }
    else
    { stop("This attribute is not a list !")}
    validObject(x)
    return(x)
  }
)
