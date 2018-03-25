/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_CategoricalParametersHandler.cpp
 *  @brief In this file we implement the parameters handlers classes for categorical
 *  mixture models
 **/

#include "../include/CategoricalModels/STK_CategoricalParametersHandler.h"

namespace STK
{

/* default constructor */
ParametersHandler<Clust::Categorical_pjk_>::ParametersHandler(int nbCluster): stat_proba_(nbCluster) {}
/* copy constructor */
ParametersHandler<Clust::Categorical_pjk_>::ParametersHandler(ParametersHandler const& model): stat_proba_(model.stat_proba_) {}
/* copy operator */
ParametersHandler<Clust::Categorical_pjk_>& ParametersHandler<Clust::Categorical_pjk_>::operator=( ParametersHandler const& other)
{ stat_proba_ = other.stat_proba_; return *this; }
/* update statistics of the parameters. */
void ParametersHandler<Clust::Categorical_pjk_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  { stat_proba_[k].update(param.proba_[k]);}
}
/* set and release the computed statistics */
void ParametersHandler<Clust::Categorical_pjk_>::setStatistics(Parameters& param)
{
  for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  {
    param.proba_[k] = stat_proba_[k].mean();
    stat_proba_[k].release();
  }
}
/* Set the computed statistics */
void ParametersHandler<Clust::Categorical_pjk_>::releaseStatistics()
{
  for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  { stat_proba_[k].release();}
}
/* Initialize the parameters of the model. */
void ParametersHandler<Clust::Categorical_pjk_>::resize(Range const& rangeModalities, Range const& rangeData )
{
  for (int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  {
    stat_proba_[k].resize(rangeModalities, rangeData);
  }
}

/* default constructor */
ParametersHandler<Clust::Categorical_pk_>::ParametersHandler(int nbCluster): stat_proba_(nbCluster) {}
/* copy constructor */
ParametersHandler<Clust::Categorical_pk_>::ParametersHandler(ParametersHandler const& model): stat_proba_(model.stat_proba_) {}
/* copy operator */
ParametersHandler<Clust::Categorical_pk_>& ParametersHandler<Clust::Categorical_pk_>::operator=( ParametersHandler const& other)
{ stat_proba_ = other.stat_proba_; return *this; }
/* update statistics of the parameters. */
void ParametersHandler<Clust::Categorical_pk_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  { stat_proba_[k].update(param.proba_[k]);}
}
/* set the computed statistics */
void ParametersHandler<Clust::Categorical_pk_>::setStatistics(Parameters& param)
{
  for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  {
    param.proba_[k] = stat_proba_[k].mean();
    stat_proba_[k].release();
  }
}
/* release the computed statistics */
void ParametersHandler<Clust::Categorical_pk_>::releaseStatistics()
{
  for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  { stat_proba_[k].release();}
}
/* Initialize the parameters of the model. */
void ParametersHandler<Clust::Categorical_pk_>::resize(Range const& rangeModalities, Range const& rangeData )
{
  for (int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
  { stat_proba_[k].resize(rangeModalities);}
}


} // namespace STK

