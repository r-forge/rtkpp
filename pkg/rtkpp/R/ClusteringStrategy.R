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
#' Create an instance of [\code{\linkS4class{ClusteringStrategy}}] class 
#'
#' A strategy is a multistage empirical process for finding a
#' good estimate in the clusering estimation process.
#' 
#' A strategy is a way to find a good estimate of the parameters of a mixture model
#' when using an EM algorithm. A ``try'' is composed of three stages 
#' \itemize{
#'   \item \code{nbInit} initializations of the \code{EM}, \code{CEM} or \code{SEM} algorithm. 
#'   \item \code{nbShortRun} short iterations of the \code{EM}, \code{CEM} or \code{SEM} algorithm.
#'   \item A long run of the \code{EM}, \code{CEM} or \code{SEM} algorithm.
#' }
#' For exemple if \code{nbInit} is 5 and \code{nbShortRun} is also 5, there will be
#' 5 packets of 5 models initialized. In each packet, the best model will be ameliorated using
#' a short run. Among the 5 models ameliorated one will be estimated until convergence
#' using a long run.
#'   
#' This process can be repeated at least \code{nbTry} times. If all the tries failed, no
#' model is returned. 
#'
#' @param nbInit Integer defining the number of initialization to try. Default value: 5.
#' @param initMethod Character string with the initialisation method.
#' @param initAlgo Character string with the algorithm to use in the initialization stage.
#' Default value: "SEM" 
#' @param nbInitIteration Integer defining the maximal number of iterations in initialization algorithm
#' if \code{initAlgo} = "EM or "CEM. This is the number of iterations if \code{initAlgo} = "SEM".
#' Default value: 20.
#' @param initEpsilon Real defining the epsilon value for the algorithm.
#' epsilon is not used by the \code{SEM} algorithm. Default value: 0.01.
#' 
#' @param nbShortRun Integer defining the number of short run to try. Default value: 5.
#' @param shortRunAlgo A character string with the algorithm to use in the short run stage
#' Default value: "EM".
#' @param nbShortIteration Integer defining the maximal number of iterations in the short runs
#' if \code{shortRunAlgo} = "EM or "CEM, or the number of iterations if \code{shortRunAlgo} = "SEM".
#' Default value: 100.
#' @param shortEpsilon Real defining the epsilon value for the algorithm.
#' epsilon is not used by the \code{SEM} algorithm. Default value: 1e-04.
#' 
#' @param longRunAlgo A character string with the algorithm to use in the long run stage
#' Default value: "EM".
#' @param nbLongIteration  Integer defining the maximal number of iterations in the short runs
#' if \code{shortRunAlgo} = "EM or "CEM, or the number of iterations if \code{shortRunAlgo} = "SEM".
#' Default value: 1000.
#' @param longEpsilon Real defining the epsilon value for the algorithm.
#' epsilon is not used by the \code{SEM} algorithm. Default value: 1e-07.
#' 
#' @param nbTry number of estimation to attempt.
#'
#' @examples
#'    clusteringStrategy()
#'    clusteringStrategy(longRunAlgo= "CEM", nbLongIteration=100)
#'    clusteringStrategy(nbTry = 3, nbInit= 1, shortRunAlgo= "SEM", nbShortIteration=100)
#'
#' @return a [\code{\linkS4class{ClusteringStrategy}}] object
#' @author Serge Iovleff
#' @export
clusteringStrategy <- function( nbTry =1
                              , nbInit= 5, initMethod="random", initAlgo= "SEM", nbInitIteration=20, initEpsilon=0.01
                              , nbShortRun= 5, shortRunAlgo= "EM", nbShortIteration=100, shortEpsilon=1e-04
                              , longRunAlgo= "EM", nbLongIteration=1000, longEpsilon=1e-07
                              )
{
  # create init
  initMethod = clusteringInit(initMethod, nbInit, initAlgo, nbInitIteration, initEpsilon);
  # create shortAlgo
  shortAlgo = clusteringAlgo(shortRunAlgo, nbShortIteration, shortEpsilon);
  # create longAlgo
  longAlgo = clusteringAlgo(longRunAlgo, nbLongIteration, longEpsilon);
  # create strategy
  new("ClusteringStrategy", nbTry =nbTry, initMethod =initMethod, shortAlgo =shortAlgo, longAlgo =longAlgo);
}

###################################################################################
#' Constructor of [\code{\linkS4class{ClusteringStrategy}}] class
#'
#' This class encapsulate the parameters of the strategy of estimation of the rtkpp
#' Clustering models.
#'
#'   @slot nbTry Integer defining the number of tries. Default value: 1.
#'   @slot initMethod A [\code{\linkS4class{ClusteringInit}}] object defining the way to
#'   initialize the estimation method. 
#'   @slot shortAlgo A [\code{\linkS4class{ClusteringAlgo}}] object defining the algorithm
#'   to use during the short runs of the estimation method. 
#'   @slot longAlgo A [\code{\linkS4class{ClusteringAlgo}}] object defining the algorithm
#'   to use during the long run of the estimation method.
#'
#' @examples
#'   new("ClusteringStrategy")
#'   new("ClusteringStrategy", shortAlgo=clusteringAlgo("SEM",1000))
#'   getSlots("ClusteringStrategy")
#'
#' @author Serge Iovleff
#' 
#' @name ClusteringStrategy
#' @rdname ClusteringStrategy-class
#' @aliases ClusteringStrategy-class
#' @exportClass ClusteringStrategy
setClass(
    Class="ClusteringStrategy",
    slots=c( nbTry = "numeric", initMethod = "ClusteringInit", shortAlgo="ClusteringAlgo", longAlgo="ClusteringAlgo" ),
    prototype=list(  nbTry = 5, initMethod = clusteringInit(), shortAlgo=clusteringAlgo("EM",100,1e-04), longAlgo=clusteringAlgo("EM",1000,1e-07)),
    # validity function
    validity=function(object)
    {
      if (round(object@nbTry)!=object@nbTry)
      {stop("nbTry must be an integer.")}
      if( object@nbTry < 1 ) # can't be zero
      {stop("nbTry must be strictly greater than 0");}
      if(class(object@initMethod)[1] != "ClusteringInit")
      {stop("shortAlgo is not of a clustering algorithm (must be an instance of the class ClusteringAlgo).")}
      if(class(object@shortAlgo)[1] != "ClusteringAlgo")
      {stop("shortAlgo is not of a clustering algorithm (must be an instance of the class ClusteringAlgo).")}
      if(class(object@longAlgo)[1] != "ClusteringAlgo")
      {stop("longAlgo is not of a clustering algorithm (must be an instance of the class ClusteringAlgo).")}
      return(TRUE)
    }
)


###################################################################################
#' Create an instance of the [\code{\linkS4class{ClusteringStrategy}}] class using new/initialize.
#' 
#' Initialization method. Used internally in the `rtkpp' package.
#' 
# @seealso \code{\link{initialize}}
#' @keywords internal
# @rdname initialize-methods
#'
setMethod(
  f="initialize",
  signature=c("ClusteringStrategy"),
  definition=function(.Object, nbTry, initMethod, shortAlgo, longAlgo)
  {
    # for nbtry
    if(missing(nbTry)) {.Object@nbTry<-5}
    else  {.Object@nbTry<-nbTry}
    # for initMethod
    if( missing(initMethod) ){ .Object@initMethod<-clusteringInit() }
    else{.Object@initMethod<-initMethod}
    # for shortAlgo
    if(missing(shortAlgo)){ .Object@shortAlgo<-clusteringAlgo("EM", 100, 1e-04) }
    else{.Object@shortAlgo<-shortAlgo}
    # for longAlgo
    if(missing(longAlgo)){ .Object@longAlgo<-clusteringAlgo("EM", 1000, 1e-07) }
    else{.Object@longAlgo<-longAlgo}
    # validate
    validObject(.Object)        
    return(.Object)
  }
)

###################################################################################
# @name print
# @docType methods
#' @aliases print print-strategy,ClusteringStrategy-method
#' @rdname print-methods
setMethod(
  f="print",
  signature=c("ClusteringStrategy"),
  function(x,...){
    cat("****************************************\n")
    cat("*** Clustering Strategy:\n")
    cat("* number of try         = ", x@nbTry, "\n")
    cat("****************************************\n")
    cat("*** Initialization :\n")
    cat("* method = ", x@initMethod@method, "\n")
    cat("* number of init       = ", x@initMethod@nbInit, "\n")
    cat("* algorithm            = ", x@initMethod@algo@algo, "\n")
    cat("* number of iterations = ", x@initMethod@algo@nbIteration, "\n")
    cat("* epsilon              = ", x@initMethod@algo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** short algorithm :\n")
    cat("* algorithm            = ", x@shortAlgo@algo, "\n")
    cat("* number of iterations = ", x@shortAlgo@nbIteration, "\n")
    cat("* epsilon              = ", x@shortAlgo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** long algorithm :\n")
    cat("* algorithm            = ", x@longAlgo@algo, "\n")
    cat("* number of iterations = ", x@longAlgo@nbIteration, "\n")
    cat("* epsilon              = ", x@longAlgo@epsilon, "\n")    
    cat("****************************************\n")
  }
)

###################################################################################
# @name show
# @docType methods
#' @rdname show-methods
#' @aliases show show-strategy,ClusteringStrategy-method
setMethod(
  f="show",
  signature=c("ClusteringStrategy"),
  function(object){
    cat("****************************************\n")
    cat("*** Clustering Strategy:\n")
    cat("* number of try         = ", object@nbTry, "\n")
    cat("****************************************\n")
    cat("*** Initialization :\n")
    cat("* method = ", object@initMethod@method, "\n")
    cat("* number of init       = ", object@initMethod@nbInit, "\n")
    cat("* algorithm            = ", object@initMethod@algo@algo, "\n")
    cat("* number of iterations = ", object@initMethod@algo@nbIteration, "\n")
    cat("* epsilon              = ", object@initMethod@algo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** short Algorithm :\n")
    cat("* algorithm            = ", object@shortAlgo@algo, "\n")
    cat("* number of iterations = ", object@shortAlgo@nbIteration, "\n")
    cat("* epsilon              = ", object@shortAlgo@epsilon, "\n")
    cat("****************************************\n")
    cat("*** long algorithm :\n")
    cat("* algorithm            = ", object@longAlgo@algo, "\n")
    cat("* number of iterations = ", object@longAlgo@nbIteration, "\n")
    cat("* epsilon              = ", object@longAlgo@epsilon, "\n")    
    cat("****************************************\n")
  }
)

###################################################################################
# @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [,ClusteringStrategy-method
#'
setMethod(
  f="[", 
  signature(x = "ClusteringStrategy"),
  definition=function(x,i,j,drop){
    if ( missing(j) ){
      switch(EXPR=i,
        "nbTry"={return(x@nbTry)},
        "initMethod"={return(x@initMethod)},
        "shortAlgo"={return(x@shortAlgo)},
        "longAlgo"={return(x@longAlgo)},
        stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
  }
)

###################################################################################
#' @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [<-,ClusteringStrategy-method
setReplaceMethod(
  f="[", 
  signature(x = "ClusteringStrategy"), 
  definition=function(x,i,j,value){
    if ( missing(j) )
    {
      switch(EXPR=i,
             "nbTry"={x@nbTry<-value},
             "initMethod"={x@initMethod<-value},
             "shortAlgo"={x@shortAlgo<-value},
             "longAlgo"={x@longAlgo<-value},
             stop("This attribute doesn't exist !")
      )
    }else{
      stop("This attribute is not a list !")
    }
    validObject(x)
    return(x)
  }
)
