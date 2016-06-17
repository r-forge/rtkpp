/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

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
 * created on: 10 May 2016
 * Author:   Iovleff, serge.iovleff@stkpp.org
 **/

/** @file ILauncher.h
 *  @brief In this file we define the interface class for launchers.
 **/


#ifndef STK_ILAUNCHER_H
#define STK_ILAUNCHER_H

#include "RDataHandler.h"

namespace STK
{

/** The ILauncher allow to create the composer or learner for estimate or
 *  learn a  mixture model with less effort.
 **/
class ILauncher: public IRunnerBase
{
  public:
    /** constructor.
     * @param model a reference on the current model
     * @param models a list of model name
     **/
    ILauncher( SEXP model, SEXP models);
    /** constructor with a list of component.
     *  @param model a reference on the current model
     **/
    ILauncher( SEXP model);
    /** destructor. */
    virtual ~ILauncher();
    /** @return the model */
    inline Rcpp::S4 const& s4_model() const { return s4_model_;}

  protected:
    /** create the data sets with real data */
    void createContinuousDataSets(std::string const& idData, std::string const& idModel,
                                  Rcpp::S4 s4_component, Clust::Mixture model);
    /** create the data sets with integer data */
    void createDiscreteDataSets(std::string const& idData, std::string const& idModel,
                                Rcpp::S4 s4_component, Clust::Mixture model);
    /** create the mixtures in the given model */
    void createMixtures(IMixtureStatModel* p_model);

    /** get the parameters */
    void getParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component);
    /** get the diagonal Gaussian parameters */
    void getDiagGaussianParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component);
    /** get the Poisson parameters */
    void getPoissonParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component);
    /** get the gamma parameters */
    void getGammaParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component);
    /** get the gamma parameters */
    void getCategoricalParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component);
    /** get the kernel parameters */
    void getKernelParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4& s4_component);

    /** diagonal Gaussian mixture models manager */
    DiagGaussianMixtureManager<RDataHandler> diagGaussianManager_;
    /** Poisson mixture models manager */
    PoissonMixtureManager<RDataHandler> poissonManager_;
    /** gamma mixture models manager */
    GammaMixtureManager<RDataHandler> gammaManager_;
    /** categorical mixture models manager */
    CategoricalMixtureManager<RDataHandler> categoricalManager_;
    /** kernel mixture models manager */
    KernelMixtureManager<RDataHandler> kernelManager_;

    /** data handler */
    RDataHandler handler_;
    /** model from the R side */
    Rcpp::S4              s4_model_;
    /** vector with the model names to try */
    Rcpp::CharacterVector v_models_;
    /** Is the model with mixed data ? */
    bool isMixedData_;
};

} // namespace STK

#endif /* STK_ILAUNCHER_H */
