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
ILauncher::ILauncher( Rcpp::S4 model, Rcpp::CharacterVector models)
                    : ILauncherBase(model)
                    , v_models_(models)
                    , handler_()
                    , diagGaussianManager_(handler_)
                    , poissonManager_(handler_)
                    , gammaManager_(handler_)
                    , categoricalManager_(handler_)
                    , kernelManager_(kerHandler_)
{}
/* facade design pattern.
 * The ILauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
ILauncher::ILauncher( Rcpp::S4 model)
                    : ILauncherBase(model)
                    , v_models_()
                    , handler_()
                    , diagGaussianManager_(handler_)
                    , poissonManager_(handler_)
                    , gammaManager_(handler_)
                    , categoricalManager_(handler_)
                    , kernelManager_(kerHandler_)
{}

/* destructor. */
ILauncher::~ILauncher()
{}


/* create the managers for models with real data */
void ILauncher::createContinuousDataSets( std::string const& idData
                                        , Rcpp::S4 s4_component
                                        , Clust::Mixture model
                                        )
{
  NumericMatrix m_data = s4_component.slot("data");
  handler_.addData(m_data, idData, Clust::mixtureToString(model));
}

/* create the managers for models with real data */
void ILauncher::createDiscreteDataSets( std::string const& idData
                                      , Rcpp::S4 s4_component
                                      , Clust::Mixture model
                                      )
{
  IntegerMatrix m_data = s4_component.slot("data");
  handler_.addData(m_data, idData, Clust::mixtureToString(model));
}

/* create the mixtures in the given learner */
void ILauncher::createMixtures(IMixtureStatModel* p_model)
{
  p_model->createMixture(diagGaussianManager_);
  p_model->createMixture(poissonManager_);
  p_model->createMixture(gammaManager_);
  p_model->createMixture(categoricalManager_);
}

/* fill the s4_component with the parameters */
void ILauncher::getParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4 s4_component)
{
  std::string rModelName = s4_component.slot("modelName");
  bool freeProp;
  switch (Clust::mixtureToMixtureClass(Clust::stringToMixture(rModelName, freeProp)))
  {
    case Clust::DiagGaussian_:
      setDiagGaussianParametersToComponent(p_model, diagGaussianManager_, idData, s4_component);
      break;
    case Clust::Poisson_:
      setPoissonParametersToComponent(p_model, poissonManager_, idData, s4_component);
      break;
    case Clust::Gamma_:
      setGammaParametersToComponent(p_model, gammaManager_, idData, s4_component);
      break;
    case Clust::Categorical_:
      setCategoricalParametersToComponent(p_model, categoricalManager_, idData, s4_component);
      break;
    case Clust::Kmm_:
      setKernelParametersToComponent(p_model, kernelManager_, idData, s4_component);
      break;
    case Clust::unknown_mixture_class_:
      break;
    default:
      break;
  }
}

} // namespace STK

