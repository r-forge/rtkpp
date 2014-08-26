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
#-----------------------------------------------------------------------
#' @include ClusterModels.R
NULL

#' Definition of the [\code{\linkS4class{ClusterDiagGaussianModel}}] class
#'
#' This class defines a diagonal Gaussian mixture Model.
#'
#' This class inherits from the [\code{\linkS4class{IClusterModel}}] class.
#' A diagonal gaussian model is a mixture model of the form:
#' \deqn{
#'   f({x}|\boldsymbol{\theta})
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \phi(x_j;\mu_{jk},\sigma^2_{jk})
#'    \quad x \in {R}^d.
#' }
#'
#' @slot meankj  Matrix with the mean of the jth variable in the kth cluster.
#' @slot sigma2kj  Matrix with the variance of the jth variable in the kth cluster.
#'
#' @examples
#' getSlots("ClusterDiagGaussianModel")
#' data(geyser)
#' new("ClusterDiagGaussianModel", data=geyser)
#'
#' @author Serge Iovleff
#'
#' @name ClusterDiagGaussianModel
#' @rdname ClusterDiagGaussianModel-class
#' @aliases ClusterDiagGaussianModel-class
#' @exportClass ClusterDiagGaussianModel
#'
setClass(
    Class="ClusterDiagGaussianModel",
    representation( meankj = "matrix", sigma2kj = "matrix"),
    contains=c("IClusterModel"),
    prototype=list( meankj   = matrix(nrow=0,ncol=0), sigma2kj = matrix(nrow=0,ncol=0) ),
    validity=function(object)
    {
      if (nrow(object@meankj)!=object@nbCluster)
      {stop("meankj must have nbCluster rows.")}
      if (ncol(object@meankj)!=ncol(object@data))
      {stop("meankj must have nbVariable columns.")}
      if (nrow(object@sigma2kj)!=object@nbCluster)
      {stop("sigma2kj must have nbCluster rows.")}
      if (ncol(object@sigma2kj)!=ncol(object@data))
      {stop("sigma2kj must have nbVariable columns.")}
      if (!(object@modelName %in% c("gaussian_sjk", "gaussian_sj", "gaussian_sk", "gaussian_s")))
      {stop("Invalid Gaussian model name.")}
      return(TRUE)
    }
)


#' Initialize an instance of the [\code{\linkS4class{ClusterDiagGaussianModel}}] class.
#'
#' Initialization method. Used internally in the `rtkpp' package.
#'
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterDiagGaussianModel"),
    definition=function(.Object, data, nbCluster=2, modelName="gaussian_sjk")
    {
      .Object <- callNextMethod(.Object, data, nbCluster, modelName)
      # resize
      nbVariable <- ncol(.Object@data)
      # for modelName
      if(missing(modelName)) {.Object@modelName<-"gaussian_sjk"}
      else  {.Object@modelName<-modelName}
      .Object@meankj <- matrix(0, nbCluster, nbVariable)
      .Object@sigma2kj <- matrix(1, nbCluster, nbVariable)
      # validate
      validObject(.Object)
      return(.Object)
    }
)

#' @rdname print-methods
#' @aliases print print,ClusterDiagGaussianModel-method
#'
setMethod(
  f="print",
  signature=c("ClusterDiagGaussianModel"),
  function(x,...){
    cat("****************************************\n")
    callNextMethod()
    cat("****************************************\n")
    for(k in 1:length(x@pk))
    {
      cat("*** Cluster: ",k,"\n")
      cat("* Proportion = ", format(x@pk[k]), "\n")
      cat("* Means      = ", format(x@meankj[k,]), "\n")
      cat("* Variances  = ", format(x@sigma2kj[k,]), "\n")
      cat("****************************************\n")
    }
  }
)

# @name show
# @docType methods
#' @rdname show-methods
#' @aliases show-diagGaussian,ClusterDiagGaussianModel,ClusterDiagGaussianModel-method
setMethod(
    f="show",
    signature=c("ClusterDiagGaussianModel"),
    function(object){
      cat("****************************************\n")
      callNextMethod()
      cat("****************************************\n")
      for(k in 1:length(object@pk))
      {
        cat("*** Cluster: ",k,"\n")
        cat("* Proportion = ", format(object@pk[k]), "\n")
        cat("* Means      = ", format(object@meankj[k,]), "\n")
        cat("* Variances  = ", format(object@sigma2kj[k,]), "\n")
        cat("****************************************\n")
      }
    }
)

#' Plotting of a class [\code{\linkS4class{ClusterDiagGaussianModel}}]
#'
#' Plotting data from a [\code{\linkS4class{ClusterDiagGaussianModel}}] object
#' using the estimated parameters and partition.
#'
#' @param x an object of class [\code{\linkS4class{ClusterDiagGaussianModel}}]
#' @param y a list of variables to plot (subset). Variables names or indices.
#' If missing all the variables are represented.
#' @param ... further arguments passed to or from other methods
#'
#' @importFrom graphics plot
#' @name plot
#' @aliases plot plot,ClusterDiagGaussianModel-method
#' @docType methods
#' @rdname plot-methods
#' @exportMethod plot
#'
#' @seealso \code{\link{plot}}
#' @examples
#'   ## for quantitative case
#'   data(iris)
#'   mixt <- clusterDiagGaussian(iris[1:4],3)
#' \dontrun{
#'   plot(mixt)
#'   plot(mixt, c(1,3))
#'   plot(mixt, c("Sepal.Length","Sepal.Width"))
#'   }

#'
setMethod(
    f="plot",
    signature=c("ClusterDiagGaussianModel"),
    function(x, y, ...)
    { # use generic plot
      .clusterPlot(x, y, .dGauss,...);
    }
)

# wrapper of dnorm
# x a vector with the point
.dGauss <- function(x, j, k, model)
{ dnorm(x, (model@meankj)[k, j] , sqrt((model@sigma2kj)[k, j]))}
