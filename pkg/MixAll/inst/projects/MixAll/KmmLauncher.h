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

/** @file KmmLauncher.h
 *  @brief In this file we define the class forKmmLauncher
 **/


#ifndef STK_KMMLAUNCHER_H
#define STK_KMMLAUNCHER_H

#include "ILauncherBase.h"

namespace STK
{

/** KmmLauncher class is an interface class allowing to create composer or
 *  learner of mixture models with less efforts.
 **/
class KmmLauncher: public ILauncherBase
{
  public:
    /** constructor for single data
     *  @param s4_model KmmModel S4 class with the result in output
     *  @param handler an instance of the KernelHandler class with the kernels to use
     *  @param nbCluster a vector with the number of clusters to test
     *  @param models a vector of string with the model names to try
     *  @param criteria selection model criteria (string)
     **/
  KmmLauncher( Rcpp::S4 s4_model
             , Rcpp::IntegerVector const& nbCluster
             , Rcpp::CharacterVector const& models
             );
//    /** constructor for mixed data
//     *  @param model ClusterMixedData S4 class
//     *  @param nbCluster a vector with the number of clusters to test
//     *  @param critName selection model criteria (string)
//     **/
//    KmmLauncher( SEXP model, SEXP nbCluster, SEXP critName );
    /** destructor. */
    virtual ~KmmLauncher();

    /** run the estimation */
    bool run();

  protected:
    /** create the kernel mixtures in the given composer. We have to
     * use a workaround for this kind of model as they need an extra parameter
     * (the dimension) to be given
     **/
    void updateMixtures(MixtureComposer* p_composer);

    /** vector with the model names to try */
    Rcpp::CharacterVector v_models_;
    /** strategy from the R side */
    Rcpp::S4              s4_strategy_;
    /** vector with the number of cluster to try */
    Rcpp::IntegerVector   v_nbCluster_;
    /** character string with the model selection criterion name */
    std::string           criterion_;

  private:
    /** Select the best model among the models and nbCluster given.
     *  @return the value of the best criteria.
     **/
    Real selectBestSingleModel();
    /** Select the best model among the models and nbCluster given.
     *  @return the value of the best criteria.
     **/
    Real selectBestMixedModel();
    /** facade with the main composer */
    MixtureComposerFacade<RDataHandler> facade_;
    /** Is the model with mixed data ? */
    bool isMixedData_;
};

} // namespace STK

#endif /* STK_IKMMLAUNCHER_H */
