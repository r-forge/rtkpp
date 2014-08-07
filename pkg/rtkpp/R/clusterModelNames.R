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
#' Build a vector of diagonal Gaussian model names.
#' 
#' In a diagonal Gaussian model, we assume that the variance
#' matrices are diagonal in each cluster. This gives rise to 8 models:
#' \enumerate{
#'  \item {The proportions can be equal or free.}
#'  \item {The variances can be equal or free for all the variables.}
#'  \item {The variances can be equal or free for all the clusters.}
#' }
#'   
#' The model names are summarized in the following array: 
#' \tabular{llll}{
#'     Model Name      \tab Prop. \tab Variance of the variables \tab Variance of the clusters \cr
#'     gaussian_p_sjk  \tab Equal \tab Free  \tab Free  \cr
#'     gaussian_p_sj   \tab Equal \tab Free  \tab Equal \cr
#'     gaussian_p_sk   \tab Equal \tab Equal \tab Free  \cr
#'     gaussian_p_s    \tab Equal \tab Equal \tab Equal \cr 
#'     gaussian_pk_sjk \tab Free  \tab Free  \tab Free  \cr
#'     gaussian_pk_sj  \tab Free  \tab Free  \tab Equal \cr
#'     gaussian_pk_sk  \tab Free  \tab Equal \tab Free  \cr
#'     gaussian_pk_s   \tab Free  \tab Equal \tab Equal \cr
#' }
#'
#' @title diagGaussianNames: get the names of the diagonal Gaussian models
#' 
#' @param prop A character string equal to "equal", "free" or "all". Default is "all".
#' @param varianceVariables A character string equal to "equal", "free" or "all". Default is "all".
#' @param varianceClusters A character string equal to "equal", "free" or "all". Default is "all".
#' 
#' @examples
#' diagGaussianNames()
#' diagGaussianNames("free", "equal", "free") # same as c("gaussian_p_sk", "gaussian_pk_sk")
#'
#' @rdname diagGaussianNames
#' @export diagGaussianNames
diagGaussianNames <- function(prop = "all", varianceVariables="all", varianceClusters = "all")
{
  if(sum(prop %in% c("equal","free","all")) != 1)
  { stop("prop is not valid. See ?diagGaussianNames for the list of prop.")}
  if(sum(varianceVariables %in% c("equal","free","all")) != 1)
  { stop("varianceVariables is not valid. See ?diagGaussianNames for the list of varianceVariables.")}
  if(sum(varianceClusters %in% c("equal","free","all")) != 1)
  { stop("varianceClusters is not valid. See ?diagGaussianNames for the list of varianceClusters.")}
  
  all = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk", "gaussian_pk_s"
         , "gaussian_p_sjk", "gaussian_p_sj", "gaussian_p_sk", "gaussian_p_s")
  propFree = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk", "gaussian_pk_s")
  propEqual = c( "gaussian_p_sjk", "gaussian_p_sj", "gaussian_p_sk", "gaussian_p_s")
  varVarFree = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_p_sjk", "gaussian_p_sj")
  varVarEqual = c( "gaussian_pk_sk", "gaussian_pk_s", "gaussian_p_sk", "gaussian_p_s")
  varClustFree = c( "gaussian_pk_sjk", "gaussian_pk_sk", "gaussian_p_sjk", "gaussian_p_sk")
  varClustEqual = c( "gaussian_pk_sj", "gaussian_pk_s", "gaussian_p_sj", "gaussian_p_s")
  
  res = all;
  if (prop == "free")  { res = intersect(res, propFree);}
  if (prop == "equal") { res = intersect(res, propEqual);}
  if (varianceVariables =="free")  { res = intersect(res, varVarFree);}
  if (varianceVariables == "equal") { res = intersect(res, varVarEqual);}
  if (varianceClusters =="free")  { res = intersect(res, varClustFree);}
  if (varianceClusters =="equal") { res = intersect(res, varClustEqual);}
  
  res
}

#' check if a vector of diagonal Gaussian model name comply 
#' @rdname checkModelNames
#' @keywords internal
checkDiagGaussianNames <- function(modelNames)
{
  nb = length(modelNames)
  if ( nb == 0 ) { return(FALSE);}
  
  all = c( "gaussian_pk_sjk", "gaussian_pk_sj", "gaussian_pk_sk", "gaussian_pk_s"
         , "gaussian_p_sjk", "gaussian_p_sj", "gaussian_p_sk", "gaussian_p_s")
  for (i in 1:nb)
  {  if ( sum(modelNames[i] %in% all) != 1 ) { return(FALSE);}
  }
  return(TRUE)
}