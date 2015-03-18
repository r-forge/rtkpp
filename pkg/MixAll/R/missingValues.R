# TODO: Add comment
#
# Author: iovleff
###############################################################################

#' Return the missing values of a component or a cluster class.
#'
#' The missing methods allow the user to get the imputed mssing
#' values from a mixture model.
#'
#' @param x an object that can return the imputed missing values
#'
#' @return A matrix with three columns (row index, column index, value)
#'
#' @name missingValues
#' @docType methods
#' @rdname missingValues-methods
#' @exportMethod missingValues
#'
#' @examples
#'   data(geyser)
#'   model <- clusterDiagGaussian(geyser,3)
#'   missingValues(model)
setGeneric(
  name = "missingValues",
  function(x)
  { standardGeneric("missingValues")}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussian-method
setMethod(
  "missingValues",
  c("ClusterHeterogeneous"),
  function(x)
  {
    nbData <- length(x@ldata)
    res <- vector("list", nbData)
    if(nbData>0)
    {
      for (l in 1:nbData)
      { res[[l]]  <- cbind(x@ldata[[l]]@missing, (x@ldata[[l]]@data)[x@ldata[[l]]@missing]);}
    }
    return(res)
  }
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussianComponent-method
setMethod(
  "missingValues",
  c("ClusterDiagGaussianComponent"),
  function(x){ return(cbind(x@missing, x@data[x@missing]));}
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussian-method
setMethod(
  "missingValues",
  c("ClusterDiagGaussian"),
  function(x){ return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterGammaComponent-method
setMethod(
    "missingValues",
    c("ClusterGammaComponent"),
    function(x){ return(cbind(x@missing, x@data[x@missing]));}
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterGamma-method
setMethod(
  "missingValues",
  c("ClusterGamma"),
  function(x){ return(missingValues(x@component));}
)


#' @rdname missingValues-methods
#' @aliases missingValues,ClusterCategoricalComponent-method
setMethod(
    f="missingValues",
    signature=c("ClusterCategoricalComponent"),
    function(x){ return(cbind(x@missing, x@data[x@missing]));}
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterCategorical-method
setMethod(
    f="missingValues",
    signature=c("ClusterCategorical"),
    function(x){ return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterPoissonComponent-method
setMethod(
    "missingValues",
    c("ClusterPoissonComponent"),
    function(x){ return(cbind(x@missing, x@data[x@missing]));}
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterPoisson-method
setMethod(
    "missingValues",
    c("ClusterPoisson"),
    function(x){ return(missingValues(x@component));}
)
