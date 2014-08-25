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
#' Create an instance of the [\code{\linkS4class{ClusterDiagGaussianModel}}] class
#'
#' This function computes an optimal diagonal Gaussian mixture model according
#' to the criteria furnished and the list of model defined in
#' [\code{modelNames}], using the
#' algorithms specified in [\code{modelNames}].
#'
#' @param data frame or matrix containing the data. Rows correspond to observations and columns correspond to variables.
#' @param nbCluster numeric listing the number of clusters.
#' @param modelNames a list of models to run. By default all diagoanl Gaussian models are estimated.
#' @param strategy a [\code{\linkS4class{ClusterStrategy}}] object containing the strategy to run.
#' clusterStrategy() method by default.
#' @param criterion character defining the criterion to select the best model. The best model is the one with
#' the lowest criterion value. Possible values: "BIC", "AIC". Default is "BIC".
#'
#' @examples
#'   ## A quantitative example with the famous geyser data set
#'   data(geyser)
#'   ## with default values
#'   clusterDiagGaussian(geyser, nbCluster=2:6)
#'   ## use graphics functions
#'   model <- clusterDiagGaussian(data=geyser, nbCluster=2:3)
#'   \dontrun{
#'   plot(model)
#'   }
#'   ## print model
#'   print(model)
#'   ## get summary
#'   summary(model)
#'
#' @return An instance of the [\code{\linkS4class{ClusterDiagGaussianModel}}] class.
#' @author Serge Iovleff
#' @export
#'
clusterDiagGaussian <- function(data, nbCluster=2, modelNames=NULL, strategy=clusterStrategy(), criterion="BIC")
{
  # check nbCluster
  nbClusterModel = length(nbCluster);
  nbClusterMin = min(nbCluster);
  nbClusterMax = max(nbCluster);
  if (nbClusterMin < 2) { stop("The number of clusters must be greater or equal to 2")}

  # check criterion
  if(sum(criterion %in% c("BIC","AIC")) != 1)
  { stop("criterion is not valid. See ?clusterDiagGaussian for the list of valid criterion")}

  # check data
  data = as.matrix(data)
  if (nrow(data) <= 3*nbClusterMax) {stop("There is not enough individuals (rows) in the data set")}
  if (ncol(data) < 1) {stop("Error: empty data set")}

  # check modelNames
  if (is.null(modelNames)) { modelNames = diagGaussianNames()}
  if (!validDiagGaussianNames(modelNames))
  { stop("modelNames is not valid. See ?diagGaussianNames for the list of valid model names")}

  # check strategy
  if(class(strategy)[1] != "ClusterStrategy")
  {stop("strategy is not a Cluster Stategy class (must be an instance of the class ClusterStrategy).")}
  validObject(strategy);

  # start estimation of the models
  model = new("ClusterDiagGaussianModel", data)
  resFlag = .Call("clusterDiagGaussianModel", model, nbCluster, modelNames, strategy, criterion, PACKAGE="rtkpp")
  # set names
  colnames(model@meankj)   <- colnames(model@data)
  colnames(model@sigma2kj) <- colnames(model@data)
  if (resFlag != 1) {cat("WARNING: An error occur during the clustering process")}
  model
}



