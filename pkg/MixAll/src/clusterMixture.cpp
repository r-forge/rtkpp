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
RcppExport SEXP clusterMixture( SEXP model, SEXP nbCluster, SEXP modelNames
                              , SEXP strategy, SEXP critName
                              , SEXP nbCore )
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
 *  @param modelNames a vector of string with the model names to try
 */
RcppExport SEXP clusterKernelMixture( SEXP model, SEXP nbCluster, SEXP modelNames
                                    , SEXP strategy, SEXP critName
                                    , SEXP nbCore )
{
  BEGIN_RCPP

#ifdef _OPENMP
  int cores = Rcpp::as<int>(nbCore);
  if (cores >= 1) { omp_set_num_threads(cores);}
#endif

  Rcpp::S4 s4_model(model);
  // build Gram matrix
  Rcpp::CharacterVector r_kernelName = s4_model.slot("kernelName");
  Rcpp::DoubleVector r_kernelParameters = s4_model.slot("kernelParameters");
  std::string kernelName = Rcpp::as<std::string>(r_kernelName[0]);
  STK::Real param1, param2;
  switch (r_kernelParameters.length())
  {
    case 0:
      param1 = 1.; param2 = 0.;
      break;
    case 1:
      param1 = r_kernelParameters[0]; param2 = 0.;
      break;
    default:
      param1 = r_kernelParameters[0]; param2 =  r_kernelParameters[1];
      break;
  }
  Rcpp::S4 s4_component = s4_model.slot("component");
  STK::RMatrix<double> data = s4_component.slot("data");
  // build gram matrix and overwrite solt data with it
  STK::Kernel::IKernelBase<STK::RMatrix<double> >* p_kernel;
  switch (STK::Kernel::stringToKernelType(kernelName))
  {
    case STK::Kernel::exponential_:
      p_kernel = new STK::Kernel::Exponential<STK::RMatrix<double> >(data, param1);
      break;
    case STK::Kernel::gaussian_:
      p_kernel = new STK::Kernel::Gaussian<STK::RMatrix<double> >(data, param1);
      break;
    case STK::Kernel::linear_:
      p_kernel = new STK::Kernel::Linear<STK::RMatrix<double> >(data);
      break;
    case STK::Kernel::polynomial_:
      p_kernel = new STK::Kernel::Polynomial<STK::RMatrix<double> >(data, param1, param2);
      break;
    case STK::Kernel::rationalQuadratic_:
      p_kernel = new STK::Kernel::RationalQuadratic<STK::RMatrix<double> >(data, param1);
      break;
    default:
      return Rcpp::wrap(false);
      break;
  }
  if (!p_kernel->run()) { delete p_kernel; return Rcpp::wrap(false);}
  s4_component.slot("data") = STK::wrap(p_kernel->gram());
  delete p_kernel;
  // create a launcher
  ClusterLauncher launcher(model, nbCluster, modelNames, strategy, critName);
  // return result
  return Rcpp::wrap(launcher.run());

  END_RCPP
}

/** @param model ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 */
RcppExport SEXP clusterMixtureHeterogene( SEXP model, SEXP nbCluster
                                        , SEXP strategy, SEXP critName, SEXP nbCore  )
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
