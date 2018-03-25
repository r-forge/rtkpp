#-----------------------------------------------------------------------
#     Copyright (C) 2012-2016  Inria
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

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussian-method
setMethod(
  f="missingValues",
  c("ClusterMixedData"),
  function(x)
  {
    nbData <- length(x@ldata)
    res <- vector("list", nbData)
    if(nbData>0)
    {
      for (l in 1:nbData)
      {
        res[[l]]  <- cbind(x@ldata[[l]]@missing, (x@ldata[[l]]@data)[x@ldata[[l]]@missing]);
        colnames(res[[l]])[3] <- "value";
      }
    }
    return(res)
  }
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussianComponent-method
setMethod(
  f="missingValues",
  c("ClusterDiagGaussianComponent"),
  function(x)
  { res = cbind(x@missing, x@data[x@missing]);
    colnames(res)[3] <- "value";
    nmiss <- nrow(x@missing)
    if (nmiss > 0) { rownames(res) <- 1:nmiss}
    return(res)
  }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterDiagGaussian-method
setMethod(
  f="missingValues",
  c("ClusterDiagGaussian"),
  function(x) { return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterGammaComponent-method
setMethod(
    f="missingValues",
    c("ClusterGammaComponent"),
    function(x)
    { res = cbind(x@missing, x@data[x@missing]);
      colnames(res)[3] <- "value";
      nmiss <- nrow(x@missing)
      if (nmiss > 0) { rownames(res) <- 1:nmiss}
      return(res)
    }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterGamma-method
setMethod(
  f="missingValues",
  c("ClusterGamma"),
  function(x) { return(missingValues(x@component));}
)


#' @rdname missingValues-methods
#' @aliases missingValues,ClusterCategoricalComponent-method
setMethod(
  f="missingValues",
  signature=c("ClusterCategoricalComponent"),
  function(x)
  { res = cbind(x@missing, x@data[x@missing]);
    colnames(res)[3] <- "value";
    nmiss <- nrow(x@missing)
    if (nmiss > 0) { rownames(res) <- 1:nmiss}
    return(res)
  }
)
#' @rdname missingValues-methods
#' @aliases missingValues,ClusterCategorical-method
setMethod(
    f="missingValues",
    signature=c("ClusterCategorical"),
    function(x) { return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterPoissonComponent-method
setMethod(
  f="missingValues",
  c("ClusterPoissonComponent"),
  function(x)
  { res = cbind(x@missing, x@data[x@missing]);
    colnames(res)[3] <- "value";
    nmiss <- nrow(x@missing)
    if (nmiss > 0) { rownames(res) <- 1:nmiss}
    return(res)
  }
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterPoisson-method
setMethod(
  f="missingValues",
  c("ClusterPoisson"),
  function(x) { return(missingValues(x@component));}
)

#' @rdname missingValues-methods
#' @aliases missingValues,ClusterPredict-method
setMethod(
  f="missingValues",
  signature=c("ClusterPredict"),
  function(x)
  { res = cbind(x@missing, x@data[x@missing]);
    colnames(res)[3] <- "value";
    nmiss <- nrow(x@missing)
    if (nmiss > 0) { rownames(res) <- 1:nmiss}
    return(res)
  }
)
#' @rdname missingValues-methods
#' @aliases missingValues,KmmComponent-method
setMethod(
    f="missingValues",
    c("KmmComponent"),
    function(x)
    {
      return(NULL)
    }
)
#' @rdname missingValues-methods
#' @aliases missingValues,KernelMixtureModel-method
setMethod(
    f="missingValues",
    c("KernelMixtureModel"),
    function(x) { return(NULL);}
)
