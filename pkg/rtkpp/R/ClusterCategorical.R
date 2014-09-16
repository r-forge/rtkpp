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
#' @include IClusterModel.R
NULL

#-----------------------------------------------------------------------
#' Definition of the [\code{\linkS4class{ClusterCategoricalModel}}] class
#'
#' This class defines a categorical mixture Model. Inherits from the
#'[\code{\linkS4class{IClusterModel}}] class. A categorical mixture model is
#' a mixture model of the form
#'
#' \deqn{
#'   f({x}|\boldsymbol{\theta}) \\
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \mathcal{M}(x_j;p_{jk},1) \\
#'    \quad {x} \in \{1,\ldots,L\}^d.
#' }
#'
#' @slot pljk  Array with the probability of the jth variable in the kth cluster
#' to be l.
#'
#' @examples
#'   getSlots("ClusterCategoricalModel")
#'   data(geyser)
#'   new("ClusterCategoricalModel", data=geyser)
#'
#' @author Serge Iovleff
#'
#' @name ClusterCategoricalModel-class
#' @rdname ClusterCategoricalModel-class
#' @aliases ClusterCategoricalModel-class
#' @exportClass ClusterCategoricalModel
#'
setClass(
    Class="ClusterCategoricalModel",
    representation( pljk = "array"),
    contains=c("IClusterModel"),
    prototype=list( pljk = array(dim=c(0,0,0))),
    validity=function(object)
    {
      dims <- dim(object@pljk)
      nbVariable <- ncol(object@data)

      if (dims[3]!=object@nbCluster)
      {stop("Third dimension in pljk must be nbCluster.")}
      if (dims[2]!=ncol(object@data))
      {stop("Second dimension in pljk must be nbVariable.")}
      if (!(object@modelName %in% c("categorical_pjk", "categorical_pk")))
      {stop("Invalid Categorical model name.")}
      return(TRUE)
    }
)

#-----------------------------------------------------------------------
#' Initialize an instance of a rtkpp class.
#'
#' Initialization method of the [\code{\linkS4class{ClusterCategoricalModel}}] class.
#' Used internally in the `rtkpp' package.
#'
#' @rdname initialize-methods
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterCategoricalModel"),
    definition=function(.Object, data, nbCluster=2, modelName="categorical_pjk")
    {
      .Object <- callNextMethod(.Object, data, nbCluster, modelName)
      # for modelName
      if(missing(modelName)) {.Object@modelName<-"categorical_pjk"}
      else  {.Object@modelName<-modelName}
      # resize
      nbVariable <- ncol(.Object@data)
      # get number of modalities
      if ( is.factor(data) ) { nbModalities <- nlevels(data)}
      else {  nbModalities <- nlevels(factor(.Object@data))}
      nbVariable <- ncol(.Object@data)
      .Object@pljk <- array(data = 1/nbModalities,dim=c(nbModalities,nbVariable,.Object@nbCluster))
      # validate
      validObject(.Object)
      return(.Object)
    }
)
