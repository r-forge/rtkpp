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

/*    Project: stkpp::MixAll
 * created on: Mar 17, 2018
 *     Author: iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
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

  for(int l=0; l<lcomponent_.length(); ++l)
  {
    Rcpp::S4 s4_component = lcomponent_[l];
    std::string idModel = s4_component.slot("modelName");
    bool freeProp;
    Clust::Mixture model           = Clust::stringToMixture(idModel, freeProp);
    Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
    std::string idData = Clust::mixtureToString(model);

    // put data set to data handler
    if ((classModel == Clust::Categorical_)||(classModel == Clust::Poisson_))
    {
      Rcpp::IntegerMatrix r_data_int = ldata_[l];
      facade_.handler_.addData(r_data_int, idData, idModel);
    }
    else
    {
      Rcpp::NumericMatrix r_data_num = ldata_[l];
      facade_.handler_.addData(r_data_num, idData, idModel);
    }
  }

  // create composer and mixtures
  int nbCluster = s4_model_.slot("nbCluster");
  facade_.p_composer_ = new MixtureComposer(nbSample, nbCluster);
  facade_.createMixtures();

  // set proportions parameters
  RVector<double> pk((SEXP)s4_model_.slot("pk"));
  facade_.setProportions(pk);

  for(int l=0; l<lcomponent_.length(); ++l)
  {
    Rcpp::S4 s4_component = lcomponent_[l];
    std::string idModel = s4_component.slot("modelName");
    bool freeProp;
    Clust::Mixture model           = Clust::stringToMixture(idModel, freeProp);
    Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
    std::string idData = Clust::mixtureToString(model);

    // get parameters from component and set them to facade_
    ArrayXX params;
    params.move(getParameters(s4_component,idData));
    if (!facade_.setParameters( idData, params)) { return false;};
  }
  // run prediction algorithm
  p_algo_->setModel(facade_.p_composer_);
  bool flag = p_algo_->run();
  // get results
  s4_clusterPredict_.slot("pk")  = Rcpp::wrap(facade_.p_composer_->pk());
  s4_clusterPredict_.slot("tik") = Rcpp::wrap(facade_.p_composer_->tik());
  s4_clusterPredict_.slot("zi")  = Rcpp::wrap(facade_.p_composer_->zi());

  Rcpp::NumericVector fi = s4_clusterPredict_.slot("lnFi");
  Rcpp::IntegerVector zi = s4_clusterPredict_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = facade_.p_composer_->computeLnLikelihood(i);
    zi[i] += (1 - baseIdx);  // set base 1 for the class labels
  }
  //
  return flag;
}


} // namespace SK
