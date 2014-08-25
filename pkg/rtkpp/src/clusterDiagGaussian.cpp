/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  rtkpp::
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file wrap.cpp
 *  @brief In this file .
 **/


#include "RTKpp.h"


using namespace Rcpp;
using namespace STK;

/** @param model ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 */
RcppExport SEXP clusterDiagGaussianModel( SEXP model, SEXP nbCluster, SEXP modelNames, SEXP strategy, SEXP criterion )
{
  BEGIN_RCPP

  // wrap S4 model and strategy with S4 Rcpp class
  S4 R_model(model);
  S4 R_strategy(strategy);
  // wrap other fields
  IntegerVector R_nbCluster(nbCluster);
  CharacterVector R_modelNames(modelNames);
  std::string R_criterion(as<std::string>(criterion));

  // wrap data matrix with Rcpp
  NumericMatrix R_data = R_model.slot("data");
  // wrap Rcpp matrix with stk++ wrapper
  RcppMatrix<double> data(R_data);
  int nbSample = R_data.rows();
  int nbVariable = R_data.cols();

  //  build handler
  RDataHandler handler;
  // free prop or not ?
  Array1D<bool> v_free;
  for (int l= 0; l <R_modelNames.size() ; ++l)
  {
    std::string idData = "model" + typeToString<int>(l);
    std::string idModel(as<std::string>(R_modelNames[l]));
    // transform R models names to stk++ models names
    // check have been done on the R side so.... Let's go
    bool free;
    Clust::Mixture model = Clust::stringToMixture(idModel, free);
    handler.addData(R_data, idData, Clust::mixtureToString(model));
    v_free.push_back(free);
  }

  // create MixtureManager
  MixtureManager<RDataHandler> manager(handler);

  // create criterion
  ICriterion* crit =0;
  if (R_criterion == "BIC") { crit = new BICCriterion();}
  if (R_criterion == "AIC") { crit = new AICCriterion();}

  // start the estimation process, should end with the best model according to the criteria
  IMixtureComposer* p_composer;
  // update model
  Real criter = Arithmetic<Real>::max();
  for (int k=0; k <R_nbCluster.length(); ++k)
  {
    int K = R_nbCluster[k];
    for (int l=0; l <R_modelNames.size(); ++l)
    {
      // create composer
      if (v_free[l])
      { p_composer = new MixtureComposer(nbSample, nbVariable, K);}
      else
      { p_composer = new MixtureComposerFixedProp(nbSample, nbVariable, K);}
      // create current mixture and register it
      std::string idData = "model" + typeToString<int>(l);
      static_cast<MixtureComposer*>(p_composer)->createMixture(manager, idData);

      // create facade and strategy
      ClusterFacade facade(p_composer);
      facade.createFullStrategy(strategy);
      // run estimation
      if (facade.run())
      {
        // convenient static cast
        MixtureComposer* p_aux = static_cast<MixtureComposer*>(p_composer);
        // compute criterion and update model
        crit->setModel(p_aux);
        if (!crit->run()) { delete p_composer; p_composer = 0; continue;}
        if (criter > crit->value())
        {
          criter = crit->value();
          R_model.slot("modelName")    = as<std::string>(R_modelNames[l]);
          R_model.slot("nbCluster")    = R_nbCluster[k];
          R_model.slot("lnLikelihood") = p_aux->lnLikelihood();
          R_model.slot("criterion")    = criter;
          R_model.slot("pk")           = wrap(p_aux->pk());
          R_model.slot("tik")          = wrap(p_aux->tik());
          R_model.slot("zi")           = wrap(p_aux->zi());
          // get parameters of the Gaussian mixture
          Array2D<Real> params;
          p_aux->getParameters(manager,idData, params);
          Array2D<Real> meankj(K, nbVariable), sigma2kj(K, nbVariable);
          for (int k=0; k<K; ++k)
          {
            meankj.row(k)   = params.row(2*k);
            sigma2kj.row(k) = params.row(2*k+1);
          }
          R_model.slot("meankj")   = wrap(meankj);
          R_model.slot("sigma2kj") = wrap(sigma2kj);
          // compute fi
          NumericVector fi = R_model.slot("fi");
          for (int i=0; i< fi.length(); ++i) { fi[i] = p_aux->likelihood(i);}
        }
      }
      // release current composer
      delete p_composer; p_composer = 0;
    }
  }
  delete crit;
  if (p_composer) delete p_composer;
  // terminate
  // return true
  return Rcpp::wrap(1);

  END_RCPP
}

