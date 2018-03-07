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

/** @file kmmMixAll.cpp
 *  @brief In this file we launch the computation for estimating a kernel mixture model.
 **/


#include "RTKpp.h"
#include "MixAll.h"
#include <MixAll/MixAll_Util.h>
#include <MixAll/KmmLauncher.h>
#include <MixAll/ClusterFacade.h>
#include <MixAll/RDataHandler.h>

/* @param model a ClusterDiagModel S4 class
 * @param nbCluster a vector with the number of clusters to test
 * @param models a vector of string with the model names to try
 */
extern "C" SEXP kmm( SEXP model, SEXP nbCluster, SEXP models, SEXP prop, SEXP nbCore )
{
  BEGIN_RCPP

#ifdef _OPENMP
  int cores = Rcpp::as<int>(nbCore);
  if (cores > 1) { omp_set_num_threads(cores);}
  else { omp_set_num_threads(1);}
#endif

  Rcpp::S4 s4_model(model);
  Rcpp::IntegerVector r_nbCluster(nbCluster);
  Rcpp::CharacterVector r_models(models);
  std::string proportion = Rcpp::as<std::string>(prop);

  Rcpp::S4 s4_component = s4_model.slot("component");

  // create kernel and eventually Gram matrix
  std::string kernelName = Rcpp::as<std::string>(s4_model.slot("kernelName"));
  Rcpp::DoubleVector r_kernelParameters  = s4_model.slot("kernelParameters");
  bool kernelComputation = Rcpp::as<bool>(s4_model.slot("kernelComputation"));

  try
  {
    STK::Kernel::IKernel* p_kernel = createKernel(s4_component, kernelName, r_kernelParameters, kernelComputation);
    if (!p_kernel) { return Rcpp::wrap(false);}

    // create kernel handler
    STK::KernelHandler handler;
    for (int i=0; i<r_models.length(); ++i)
    {
      std::string idData  = "idData" + STK::typeToString(i);
      std::string idModel  = Rcpp::as<std::string>(r_models[i]);
      handler.addKernel(p_kernel, idData, idModel);
    }
    // create launcher
    STK::KmmLauncher launcher(s4_model, handler, r_nbCluster, r_models);
    // return result
    //return Rcpp::wrap(launcher.run());
    return Rcpp::wrap(true);
  }
  catch (STK::Exception const& e)
  {
    return Rcpp::wrap(false);
  };

  END_RCPP
}

