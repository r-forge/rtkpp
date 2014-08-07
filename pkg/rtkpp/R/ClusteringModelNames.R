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

###################################################################################
#' Interface base class [\code{\linkS4class{IClusteringModelNames}}]
#'
#' This class is a virtual class of the Clustering model names.
#'
#' 
#' @slot listModels character list of models.
#' @slot typeModels character string with the type of model selected. Possible values:
#' "equal", "free", all". Default value: "all".
#'
#' @examples
#'   getSlots("IClusteringModelNames")
#'
#' @author Serge Iovleff
#' 
#' @name IClusteringModelNames
#' @rdname IClusteringModelNames-class
#' @aliases IClusteringModelNames-class
#' @exportClass IClusteringModelNames
#'
setClass(
    Class="IClusteringModelNames",
    slots=c(listModels = "character", typeModels = "character"),
    prototype=list(listModels = character(1), typeModels = "all"),
    # validity function
    validity=function(object)
    {
      if(sum(object@typeModels %in% c("equal","free","all")) != 1)
      {  stop("typeModels is not valid. See ?clusteringModelNames for the list of available models.")}
      return(TRUE)
    }
)

###################################################################################
# @name print
# @docType methods
#' @aliases print print-modelnames,IClusteringModelNames-method
#' @rdname print-methods
setMethod(
  f="print",
  signature=c("IClusteringModelNames"),
  function(x,...){
    cat("****************************************\n")
    cat("*** Models:\n")
    cat("* list = ", x@listModels, "\n")
    if ( x@typeModels == "all" )
      cat("* This list includes models with free and equal proportions.\n")
    else if ( x@typeModels == "free" )
      cat("* This list includes only models with free proportions.\n")
    else if ( x@typeModels == "equal" )
      cat("* This list includes only models with equal proportions.\n")
    cat("****************************************\n")
  }
)

###################################################################################
# @name show
# @docType methods
#' @rdname show-methods
#' @aliases show show-modelnames,IClusteringModelNames-method
setMethod(
  f="show",
  signature=c("IClusteringModelNames"),
  function(object){
    cat("****************************************\n")
    cat("*** Models:\n")
    cat("* list = ", object@listModels, "\n")
    if ( object@typeModels == "all" )
      cat("* This list includes models with free and equal proportions.\n")
    else if ( object@typeModels == "free" )
      cat("* This list includes only models with free proportions.\n")
    else if ( object@typeModels == "equal" )
      cat("* This list includes only models with equal proportions.\n")
    cat("****************************************\n")
  }
)

