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

#' Create a vector of Kernel mixture model (KMM) names.
#'
#' In a Kernel mixture model, we can assume that
#' \enumerate{
#'  \item {The standard deviations are equal or free for all the clusters.}
#' }
#' This give rise to two models.
#'
#' The model names are summarized in the following array:
#' \tabular{ll}{
#'  Model Name  \tab standard deviation between clusters \cr
#'  kmm_sk      \tab Free                                \cr
#'  kmm_s       \tab Equal                               \cr
#' }
#'
#' @param sdBetweenCluster A character string equal to "equal", "free" or "all". Default is "all".
#'
#' @return A vector of character with the model names.
#' @examples
#' kmmNames()
#' ## same as c("kmm_sk")
#' kmmNames( sdBetweenCluster= "free")
#'
#' @rdname kmmNames
#' @export
kmmNames <- function( sdBetweenCluster = "all")
{
  if(sum(sdBetweenCluster %in% c("equal","free","all")) != 1)
  { stop("sdBetweenCluster is not valid. See ?kmmNames for the list of sdBetweenCluster.")}

  all = c( "kmm_sk", "kmm_s")
  sdFree   = c( "kmm_sk")
  sdEqual  = c( "kmm_s")

  res = all;
  if (sdBetweenCluster =="free")  { res = intersect(res, sdFree);}
  if (sdBetweenCluster =="equal") { res = intersect(res, sdEqual);}

  res
}

#' check if a vector of kernel mixture model (KMM) name is correct.
#' @param names a vector of character with the names to check
#' @return TRUE if the names in the vector names are valid, FALSE otherwise.
#' @rdname kmmNames
kmmValidModelNames <- function(names)
{
  nb = length(names)
  if ( nb == 0 ) { return(FALSE);}

  all = kmmNames()
  for (i in 1:nb)
  { if ( sum(names[i] %in% all) != 1 ) { return(FALSE);}}
  return(TRUE)
}

#' check if a vector of kernel name is correct.
#' @rdname kmmNames
kmmValidKernelNames <- function(names)
{
  nb = length(names)
  if ( nb == 0 ) { return(FALSE);}
    
  all = c("GAUSSIAN","POLYNOMIAL", "LAPLACE","LINEAR","RATIONALQUADRATIC","HAMMING")
  for (i in 1:nb)
  { 
    if ( sum(toupper(names[i]) %in% all) != 1 ) { return(FALSE);}
  }
  return(TRUE)
}

