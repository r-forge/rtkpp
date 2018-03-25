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

/** @file STK_PoissonParametersHandler.cpp
 *  @brief In this file we implement the parameters handlers classes for Poisson
 *  mixture models
 **/

#include "../include/PoissonModels/STK_PoissonParametersHandler.h"

namespace STK
{

/* default constructor. All lambdas are initialized to 1. */
ParametersHandler<Clust::Poisson_ljlk_>::ParametersHandler( int nbCluster)
                        : stat_lambdak_(nbCluster), stat_lambdaj_()
{}
/* copy constructor.
 * @param param the parameters to copy.
 **/
ParametersHandler<Clust::Poisson_ljlk_>::ParametersHandler( ParametersHandler const& param)
                        : stat_lambdak_(param.stat_lambdak_)
                        , stat_lambdaj_(param.stat_lambdaj_)
{}
/* destructor */
ParametersHandler<Clust::Poisson_ljlk_>::~ParametersHandler() {}
/* copy operator */
ParametersHandler<Clust::Poisson_ljlk_>& ParametersHandler<Clust::Poisson_ljlk_>::operator=( ParametersHandler const& other)
{
  stat_lambdaj_ = other.stat_lambdaj_;
  stat_lambdak_ = other.stat_lambdak_;
  return *this;
}
/* update statistics of the parameters. */
void ParametersHandler<Clust::Poisson_ljlk_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
  { stat_lambdak_[k].update(param.lambdak_[k]);}
  stat_lambdaj_.update(param.lambdaj_);
}
/* Set the computed statistics */
void ParametersHandler<Clust::Poisson_ljlk_>::setStatistics(Parameters& param)
{
  for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
  {
    param.lambdak_[k] = stat_lambdak_[k].mean();
    stat_lambdak_[k].release();
  }
  param.lambdaj_ = stat_lambdaj_.mean();
  stat_lambdaj_.release();

}
/* Release the computed statistics */
void ParametersHandler<Clust::Poisson_ljlk_>::releaseStatistics()
{
  for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
  { stat_lambdak_[k].release();}
  stat_lambdaj_.release();
}
/* Initialize the statistics. */
void ParametersHandler<Clust::Poisson_ljlk_>::resize(Range const& range)
{
  stat_lambdaj_.resize(range);
  for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
  { stat_lambdak_[k].release();}
}


/* default constructor. All lambdas are initialized to 1. */
ParametersHandler<Clust::Poisson_ljk_>::ParametersHandler( int nbCluster)
                        : stat_lambda_(nbCluster)
{}
/* copy constructor.
 * @param param the parameters to copy.
 **/
ParametersHandler<Clust::Poisson_ljk_>::ParametersHandler( ParametersHandler const& param)
                        : stat_lambda_(param.stat_lambda_)
{}
/* destructor */
ParametersHandler<Clust::Poisson_ljk_>::~ParametersHandler() {}
/* copy operator */
ParametersHandler<Clust::Poisson_ljk_>& ParametersHandler<Clust::Poisson_ljk_>::operator=( ParametersHandler const& other)
{
  stat_lambda_ = other.stat_lambda_;
  return *this;
}
/* update statistics of the parameters. */
void ParametersHandler<Clust::Poisson_ljk_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].update(param.lambda_[k]);}
}
/* Set the computed statistics */
void ParametersHandler<Clust::Poisson_ljk_>::setStatistics(Parameters& param)
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  {
    param.lambda_[k] = stat_lambda_[k].mean();
    stat_lambda_[k].release();
  }
}
/* Release the computed statistics */
void ParametersHandler<Clust::Poisson_ljk_>::releaseStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].release();}
}
/* Initialize the statistics. */
void ParametersHandler<Clust::Poisson_ljk_>::resize(Range const& range)
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].resize(range);}
}

/* default constructor. All lambdas are initialized to 1. */
ParametersHandler<Clust::Poisson_lk_>::ParametersHandler( int nbCluster)
                        : stat_lambda_(nbCluster)
{}
/* copy constructor.
 * @param param the parameters to copy.
 **/
ParametersHandler<Clust::Poisson_lk_>::ParametersHandler( ParametersHandler const& param)
                        : stat_lambda_(param.stat_lambda_)
{}
/* destructor */
ParametersHandler<Clust::Poisson_lk_>::~ParametersHandler() {}
/* copy operator */
ParametersHandler<Clust::Poisson_lk_>& ParametersHandler<Clust::Poisson_lk_>::operator=( ParametersHandler const& other)
{
  stat_lambda_ = other.stat_lambda_;
  return *this;
}
/* update statistics of the parameters. */
void ParametersHandler<Clust::Poisson_lk_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].update(param.lambda_[k]);}
}
/* Set the computed statistics */
void ParametersHandler<Clust::Poisson_lk_>::setStatistics(Parameters& param)
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  {
    param.lambda_[k] = stat_lambda_[k].mean();
    stat_lambda_[k].release();
  }
}
/* Release the computed statistics */
void ParametersHandler<Clust::Poisson_lk_>::releaseStatistics()
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].release();}
}
/* Initialize the statistics. */
void ParametersHandler<Clust::Poisson_lk_>::resize(Range const& range)
{
  for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
  { stat_lambda_[k].release();}
}


} // namespace STK

