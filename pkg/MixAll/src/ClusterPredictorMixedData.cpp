/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2018  Serge Iovleff, Université Lille 1, Inria

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

/* Project: stkpp::MixAll
 * created on: Mar 17, 2018
 * Author: iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file ClusterPredictorMixedData.cpp
 *  @brief In this file we implement the ClusterPredictorMixedData class
 **/


#include "../inst/projects/MixAll/ClusterPredictorMixedData.h"

namespace STK
{
ClusterPredictorMixedData::ClusterPredictorMixedData( Rcpp::S4 s4_model, Rcpp::S4 s4_clusterPredict)
                                                    : IClusterPredictor(s4_model, s4_clusterPredict)
                                                    , lcomponent_(s4_model_.slot("lcomponent"))
                                                    , ldata_(s4_clusterPredict_.slot("ldata"))
{}
ClusterPredictorMixedData::~ClusterPredictorMixedData() {}

/* run the estimation */
bool ClusterPredictorMixedData::run()
{
  int nbSample = s4_clusterPredict_.slot("nbSample");
  // put data set to data handler
  for(int l=0; l<lcomponent_.length(); ++l)
  {
    Rcpp::S4 s4_component = lcomponent_[l];
    String idModel = s4_component.slot("modelName");
    bool freeProp;
    Clust::Mixture model           = Clust::stringToMixture(idModel, freeProp);
    Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
    String idData = Clust::mixtureToString(model);

#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In ClusterPredictorMixedData::run. Set data for idModel =") << idModel << _T("\n");
    stk_cout << _T("In ClusterPredictorMixedData::run. Set data for idData =")  << idData << _T("\n");
#endif
    if ((classModel == Clust::Categorical_)||(classModel == Clust::Poisson_))
    {
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("ClusterPredictorMixedData::run. Adding r_data_int to data handler\n");
#endif
      Rcpp::IntegerMatrix r_data_int = ldata_[l];
      facade_.handler().addData(r_data_int, idData, idModel);
    }
    else
    {
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("adding r_data_num\n");
#endif
      Rcpp::NumericMatrix r_data_num = ldata_[l];
      facade_.handler().addData(r_data_num, idData, idModel);
    }
  }

  // create composer and mixtures
  int nbCluster = s4_model_.slot("nbCluster");
  facade_.createComposer(nbSample, nbCluster);
  facade_.createMixtures();

  // set proportions parameters
  RVector<double> pk((SEXP)s4_model_.slot("pk"));
  facade_.setProportions(pk);

  for(int l=0; l<lcomponent_.length(); ++l)
  {
    Rcpp::S4 s4_component = lcomponent_[l];
    String idModel = s4_component.slot("modelName");
    String idData = Clust::mixtureToString(Clust::stringToMixture(idModel));

#ifdef STK_MIXTURE_VERBOSE
    stk_cout << _T("In ClusterPredictorMixedData::run. Set Parameters for idModel =") << idModel << _T("\n");
    stk_cout << _T("In ClusterPredictorMixedData::run. Set Parameters for idData =")  << idData << _T("\n");
#endif
    // get parameters from component and set them to facade_
    ArrayXX params;
    params.move(getParameters(s4_component,idData));
#ifdef STK_MIXTURE_VERY_VERBOSE
      stk_cout << _T("params =")  << params;
#endif
    if (!facade_.setParameters( idData, params)) { return false;};

  }
  // run prediction algorithm
  p_algo_->setModel(facade_.p_composer());
  bool flag = p_algo_->run();
  // get results
  s4_clusterPredict_.slot("pk")  = Rcpp::wrap(facade_.p_composer()->pk());
  s4_clusterPredict_.slot("tik") = Rcpp::wrap(facade_.p_composer()->tik());
  s4_clusterPredict_.slot("zi")  = Rcpp::wrap(facade_.p_composer()->zi());

  Rcpp::NumericVector fi = s4_clusterPredict_.slot("lnFi");
  Rcpp::IntegerVector zi = s4_clusterPredict_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = facade_.p_composer()->computeLnLikelihood(i);
    zi[i] += (1 - baseIdx);  // set base 1 for the class labels
  }
  // get data set from data handler
  s4_clusterPredict_.slot("ldata") = Rcpp::clone(facade_.handler().data());
  //
  return flag;
}


} // namespace SK
