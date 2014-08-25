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
#' Interface Class [\code{\linkS4class{IClusterModel}}] for Cluster models.
#'
#' This class encapsulate the common parameters of all the Cluster models.
#'
#' A Cluster model is a model of the form
#' \deqn{
#'   f({x}|\boldsymbol{\theta}) \\
#'     \sum_{k=1}^K p_k f({x};\boldsymbol{\lambda}_k,\boldsymbol{\alpha}) \\
#'    \quad {x} \in J.
#' }
#'
#' @slot data \code{\link{matrix}} with the data set to cluster.
#' @slot nbCluster Integer with the number of cluster of the model.
#' @slot pk        Vector with the proportions of each mixture.
#' @slot tik       Matrix with the posterior probability of the ith individual to belong to kth cluster.
#' @slot fi        Vector with the mixture probabilities of the individuals.
#' @slot zi        Vector of integer with the class label of the individuals.
#' @slot lnLikelihood Real given the ln-liklihood of the Cluster model.
#' @slot criterion Real given the value of the AIC or BIC criterion.
#' @slot modelName mixture model name.
#'
#' @examples
#'   getSlots("IClusterModel")
#'
#' @author Serge Iovleff
#'
#' @name IClusterModel
#' @rdname IClusterModel-class
#' @aliases IClusterModel-class
#' @exportClass IClusterModel
setClass(
  Class="IClusterModel",
  representation( data = "matrix", nbCluster = "numeric"
                , pk = "numeric", tik = "matrix"
                , fi = "numeric", zi = "integer"
                , lnLikelihood = "numeric", criterion = "numeric"
                , modelName = "character"
                , "VIRTUAL"
                ),
  prototype=list( data = matrix(nrow=0,ncol=0)
                , nbCluster = 0
                , pk = vector("numeric")
                , tik = matrix(nrow=0, ncol=0)
                , fi  = vector("numeric")
                , zi  = vector("integer")
                , lnLikelihood = -Inf
                , criterion = -Inf
                , modelName = character(1)
  ),
  # validity function
  validity=function(object)
  {
    nbSample = nrow(object@data)
    # check nbCluster dimensions
    if (round(object@nbCluster)!=object@nbCluster)
    {stop("nbCluster must be an integer.")}
    if( object@nbCluster < 2 )
    {  stop("nbCluster must be greater than 1.")}
    if (length(object@pk)!=object@nbCluster)
    {stop("pk must have length nbCluster.")}
    if (ncol(object@tik)!=object@nbCluster)
    {stop("tik must have nbCluster columns.")}
    if (nrow(object@tik)!=nbSample)
    {stop("tik must have nbSample rows.")}
    if (length(object@fi)!=nbSample)
    {stop("fi must have nbSample size.")}
    if (length(object@zi)!=nbSample)
    {stop("zi must have nbSample size.")}
    return(TRUE)
  }
)

###################################################################################
#' Initialize an instance of the [\code{\linkS4class{IClusterModel}}] class.
#'
#' Initialization method. Used internally in the `rtkpp' package.
#'
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("IClusterModel"),
    definition=function(.Object, data, nbCluster, modelName)
    {
      # for data
      if(missing(data)) {stop("data is mandatory.")}
      .Object@data<-as.matrix(data)
      # for nbCluster
      if(missing(nbCluster)) {stop("nbCluster is mandatory.")}
      .Object@nbCluster<-nbCluster
      # for nbCluster
      if(missing(modelName)) {stop("modelName is mandatory.")}
      .Object@modelName<-modelName
      # resize
      nbSample    <- nrow(.Object@data)
      nbVariable  <- ncol(.Object@data)
      .Object@pk  <- rep(1/.Object@nbCluster, nbCluster)
      .Object@tik <- matrix(1/nbCluster, nbSample, nbCluster)
      .Object@fi  <- rep(0, nbSample)
      .Object@zi  <- as.integer(rep(1, nbSample))
      # validate
      # validObject(.Object) will be called at the end of the intialization process
      return(.Object)
    }
)


###################################################################################
#' Definition of the [\code{\linkS4class{ClusterGammaModel}}] class
#'
#' This class defines a gamma mixture Model. Inherits from the
#'[\code{\linkS4class{IClusterModel}}] class. A gamma mixture model is
#' a mixture model of the form
#'
#' \deqn{
#'   f({x}|\boldsymbol{\theta}) \\
#'   =\sum_{k=1}^K p_k \prod_{j=1}^d \gamma(x_j;a_{jk},b_{jk}) \\
#'    \quad {x} \in {R}^d.
#' }
#'
#' @slot ajk  Matrix with the shape of the jth variable in the kth cluster.
#' @slot bjk  Matrix with the scale of the jth variable in the kth cluster.
#'
#' @examples
#'   getSlots("ClusterGammaModel")
#'   data(geyser)
#'   new("ClusterGammaModel", data=geyser)
#'
#' @author Serge Iovleff
#'
#' @name ClusterGammaModel
#' @rdname ClusterGammaModel-class
#' @aliases ClusterGammaModel-class
#' @exportClass ClusterGammaModel
#'
setClass(
    Class="ClusterGammaModel",
    representation( ajk = "matrix", bjk = "matrix"),
    contains=c("IClusterModel"),
    prototype=list( ajk = matrix(nrow=0,ncol=0)
                  , bjk = matrix(nrow=0,ncol=0)
                  ),
    validity=function(object)
    {
      if (ncol(object@ajk)!=object@nbCluster)
      {stop("ajk must have nbCluster columns.")}
      if (nrow(object@ajk)!=ncol(object@data))
      {stop("ajk must have nbVariable rows.")}
      if (ncol(object@bjk)!=object@nbCluster)
      {stop("bjk must have nbCluster columns.")}
      if (nrow(object@bjk)!=ncol(object@data))
      {stop("bjk must have nbVariable rows.")}
      if (!(object@modelName %in% c("gamma_ajk_bjk", "gamma_ajk_bj")))
      {stop("Invalid Gamma model name.")}
      return(TRUE)
    }
)

###################################################################################
#' Initialize an instance of the [\code{\linkS4class{ClusterGammaModel}}] class.
#'
#' Initialization method. Used internally in the `rtkpp' package.
#'
#' @keywords internal
setMethod(
    f="initialize",
    signature=c("ClusterGammaModel"),
    definition=function(.Object, data, nbCluster=2, modelName="gamma_ajk_bjk")
    {
      .Object <- callNextMethod(.Object, data, nbCluster, modelName)
      # for modelName
      if(missing(modelName)) {.Object@modelName<-"gamma_ajk_bjk"}
      else  {.Object@modelName<-modelName}
      # resize
      nbVariable <- ncol(.Object@data)
      .Object@ajk <- matrix(0, nbVariable, nbCluster)
      .Object@bjk <- matrix(1, nbVariable, nbCluster)
      # validate
      validObject(.Object)
      return(.Object)
    }
)

###################################################################################
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

###################################################################################
#' Initialize an instance of the [\code{\linkS4class{ClusterCategoricalModel}}] class.
#'
#' Initialization method. Used internally in the `rtkpp' package.
#'
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

###################################################################################
# @name [
# @docType methods
#' @rdname extract-methods
#' @aliases [,IClusterModel-method
#'
setMethod(
    f="[",
    signature(x = "IClusterModel"),
    definition=function(x,i,j,drop){
      if ( missing(j) ){
        switch(EXPR=i,
            "data"={return(x@data)},
            "nbCluster"={return(x@nbCluster)},
            "pk"={return(x@pk)},
            "tik"={return(x@tik)},
            "fi"={return(x@fi)},
            "zi"={return(x@zi)},
            "lnLikelihood"={return(x@lnLikelihood)},
            "criterion"={return(x@criterion)},
            "modelName"={return(x@modelName)},
            stop("This attribute doesn't exist !")
        )
      }else{
        stop("This attribute is not a list !")
      }
    }
)

###################################################################################
#' @rdname print-methods
#' @aliases print print,IClusterModel-method
#'
setMethod(
  f="print",
  signature=c("IClusterModel"),
  function(x,...)
  {
    cat("* nbCluster    = ", x@nbCluster, "\n")
    cat("* lnLikelihood = ", x@lnLikelihood,"\n")
    cat("* criterion    = ", x@criterion, "\n")
    cat("* model name   = ", x@modelName, "\n")
  }
)

###################################################################################
#' @rdname print-methods
#' @aliases print print,IClusterModel-method
#'
setMethod(
    f="show",
    signature=c("IClusterModel"),
    function(object)
    {
      cat("* nbCluster    = ", object@nbCluster, "\n")
      cat("* lnLikelihood = ", object@lnLikelihood,"\n")
      cat("* criterion    = ", object@criterion, "\n")
      cat("* model name   = ", object@modelName, "\n")
    }
)
