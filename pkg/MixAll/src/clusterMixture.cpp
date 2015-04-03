/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file clusterMixture.cpp
 *  @brief In this file we launch the computation for estimating a mixture model.
 **/


#include "RTKpp.h"
#include "ClusterFacade.h"
#include "ClusterLauncher.h"
#include "RDataHandler.h"

/** @param model ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param modelNames a vector of string with the model names to try
 */
RcppExport SEXP clusterMixture( SEXP model, SEXP nbCluster, SEXP modelNames, SEXP strategy, SEXP critName, SEXP nbCore )
{
  BEGIN_RCPP

#ifdef _OPENMP
  int cores = Rcpp::as<int>(nbCore);
  if (cores >= 1) { omp_set_num_threads(cores);}
#endif

#if defined(__sun) && defined(__SVR4) && defined(__sparc)
  Rcpp::CharacterVector v_modelNames(modelNames);
  std::string idModel(Rcpp::as<std::string>(v_modelNames[0]));
  bool freeProp;
  STK::Clust::Mixture mixtModel = STK::Clust::stringToMixture(idModel, freeProp);
  STK::Clust::MixtureClass modelClass = STK::Clust::mixtureToMixtureClass(mixtModel);
  // Poisson models not working with sparc and solaris
  if (modelClass == STK::Clust::Poisson_) return Rcpp::wrap(false);
#endif
  // create a launcher
  ClusterLauncher launcher(model, nbCluster, modelNames, strategy, critName);
  // return result
  return Rcpp::wrap(launcher.run());

  END_RCPP
}

/** @param model ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 */
RcppExport SEXP clusterMixtureHeterogene( SEXP model, SEXP nbCluster, SEXP strategy, SEXP critName, SEXP nbCore  )
{
  BEGIN_RCPP

#ifdef _OPENMP
  int cores = Rcpp::as<int>(nbCore);
  if (cores >= 1) { omp_set_num_threads(cores);}
#endif
  // create a launcher
  ClusterLauncher launcher(model, nbCluster, strategy, critName);
  // return result
  return Rcpp::wrap(launcher.run());

  END_RCPP
}
