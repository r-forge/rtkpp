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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 15 may 2016
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_ILauncher.cpp
 *  @brief In this file we implement the ILauncher which
 *  construct properly a mixture model.
 **/


#include "../inst/projects/MixAll/ILauncher.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The ILauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ILauncher::ILauncher( SEXP model, SEXP models)
                            : IRunnerBase()
                            , s4_model_(model)
                            , v_models_(models)
                            , handler_()
                            , diagGaussianManager_(handler_)
                            , poissonManager_(handler_)
                            , gammaManager_(handler_)
                            , categoricalManager_(handler_)
                            , kernelManager_(handler_)
                            , isMixedData_(false)
{}
/* facade design pattern.
 * The ILauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ILauncher::ILauncher( SEXP model)
                            : IRunnerBase()
                            , s4_model_(model)
                            , v_models_()
                            , handler_()
                            , diagGaussianManager_(handler_)
                            , poissonManager_(handler_)
                            , gammaManager_(handler_)
                            , categoricalManager_(handler_)
                            , kernelManager_(handler_)
                            , isMixedData_(true)
{}
/* destructor. */
ILauncher::~ILauncher()
{
}


/* create the managers for models with real data */
void ILauncher::createContinuousDataSets( std::string const& idData
                                              , std::string const& idModel
                                              , Rcpp::S4 s4_component
                                              , Clust::Mixture model
                                              )
{
  NumericMatrix m_data = s4_component.slot("data");
  RMatrix<double> data(m_data);
  handler_.addData(m_data, idData, Clust::mixtureToString(model));
}
/* create the managers for models with real data */
void ILauncher::createDiscreteDataSets( std::string const& idData
                                          , std::string const& idModel
                                          , Rcpp::S4 s4_component
                                          , Clust::Mixture model
                                          )
{
  IntegerMatrix m_data = s4_component.slot("data");
  RMatrix<int> data(m_data);
  handler_.addData(m_data, idData, Clust::mixtureToString(model));
}

/* create the mixtures in the given learner */
void ILauncher::createMixtures(IMixtureStatModel* p_model)
{
  p_model->createMixture(diagGaussianManager_);
  p_model->createMixture(poissonManager_);
  p_model->createMixture(gammaManager_);
  p_model->createMixture(categoricalManager_);
  p_model->createMixture(kernelManager_);
}

/* fill the s4_component with the parameters */
void ILauncher::getParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component)
{
  std::string rModelName = s4_component.slot("modelName");
  bool freeProp;
  switch (Clust::mixtureToMixtureClass(Clust::stringToMixture(rModelName, freeProp)))
  {
    case Clust::Gaussian_:
      getDiagGaussianParameters(p_model, idData, s4_component);
      break;
    case Clust::Poisson_:
      getPoissonParameters(p_model, idData, s4_component);
      break;
    case Clust::Gamma_:
      getGammaParameters(p_model, idData, s4_component);
      break;
    case Clust::Categorical_:
      getCategoricalParameters(p_model, idData, s4_component);
      break;
    case Clust::Kernel_:
      getKernelParameters(p_model, idData, s4_component);
      break;
    case Clust::unknown_mixture_class_:

    default:
      break;
  }
}

/* get the diagonal Gaussian parameters */
void ILauncher::getDiagGaussianParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(diagGaussianManager_,idData, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  ArrayXX mean(K, nbVariable), sigma(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    mean.row(k)  = params.row(2*k);
    sigma.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_component.slot("mean")  = Rcpp::wrap(mean);
  s4_component.slot("sigma") = Rcpp::wrap(sigma);
  // get data
  s4_component.slot("data") = diagGaussianManager_.getData<double>(idData).matrix();
}

/* get the kernel parameters */
void ILauncher::getKernelParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component)
{
  // get parameters
  ArrayXX param;
  p_model->getParameters(kernelManager_,idData, param);
  // save results in s4_model
  s4_component.slot("sigma") = Rcpp::wrap(param.col(0));
  s4_component.slot("dim")   = Rcpp::wrap(param.col(1));
  // get data -- not necessary for kernels--
  //s4_component.slot("data") = kernelManager_.getData<double>(idData).matrix();
}

/* get the diagonal Gaussian parameters */
void ILauncher::getPoissonParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(poissonManager_,idData, params);
  // save results in s4_model
  s4_component.slot("lambda")  = Rcpp::wrap(params);
  // get data
  s4_component.slot("data") = poissonManager_.getData<double>(idData).matrix();
}

/* get the gamma parameters */
void ILauncher::getGammaParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(gammaManager_,idData, params);
  // get dimensions
  int K = params.sizeRows()/2, nbVariable = params.sizeCols();
  // get results
  ArrayXX shape(K, nbVariable), scale(K, nbVariable);
  for (int k=0; k<K; ++k)
  {
    shape.row(k) = params.row(2*k);
    scale.row(k) = params.row(2*k+1);
  }
  // save results in s4_model
  s4_component.slot("shape") = Rcpp::wrap(shape);
  s4_component.slot("scale") = Rcpp::wrap(scale);
  // get data
  s4_component.slot("data") = gammaManager_.getData<double>(idData).matrix();
}

/* get the diagonal Categorical parameters */
void ILauncher::getCategoricalParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component)
{
  // get parameters
  ArrayXX params;
  p_model->getParameters(categoricalManager_,idData, params);
  params.shift(0,0);
  // save results in s4_model
  s4_component.slot("plkj") = Rcpp::wrap(params);
  // get data
  s4_component.slot("data") = categoricalManager_.getData<int>(idData).matrix();
}

} // namespace STK

