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

/** @file STK_DiagGaussianParametersHandler.cpp
 *  @brief In this file we implement the parameters handler classes for DiagGaussian
 *  mixture models
 **/

#include "../include/DiagGaussianModels/STK_DiagGaussianParametersHandler.h"

namespace STK
{
/* default constructor. */
ParametersHandler<Clust::Gaussian_sjk_>::ParametersHandler( int nbCluster)
                                                          : stat_mean_(nbCluster)
                                                          , stat_sigma_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/

ParametersHandler<Clust::Gaussian_sjk_>::ParametersHandler( ParametersHandler const& param)
                        : stat_mean_(param.stat_mean_)
                        , stat_sigma_(param.stat_sigma_)
{}
/* destructor */

ParametersHandler<Clust::Gaussian_sjk_>::~ParametersHandler() {}
/** copy operator */

ParametersHandler<Clust::Gaussian_sjk_>& ParametersHandler<Clust::Gaussian_sjk_>::operator=( ParametersHandler const& other)
{
  stat_mean_ = other.stat_mean_;
  stat_sigma_ = other.stat_sigma_;
  return *this;
}
/* update statistics of the parameters. */

void ParametersHandler<Clust::Gaussian_sjk_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].update(param.mean_[k]);
    stat_sigma_[k].update(param.sigma_[k]);
  }
}
/* Set the computed statistics */
void ParametersHandler<Clust::Gaussian_sjk_>::setStatistics(Parameters& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    param.mean_[k] = stat_mean_[k].mean();
    stat_mean_[k].release();
    param.sigma_[k] = stat_sigma_[k].mean();
    stat_sigma_[k].release();
  }
}
/* Release the computed statistics */
void ParametersHandler<Clust::Gaussian_sjk_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].release();
    stat_sigma_[k].release();
  }
}
/* Initialize the statistics. */
void ParametersHandler<Clust::Gaussian_sjk_>::resize(Range const& range)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].resize(range);
    stat_sigma_[k].resize(range);
  }
}


/* default constructor. */
ParametersHandler<Clust::Gaussian_sk_>::ParametersHandler( int nbCluster)
                        : stat_mean_(nbCluster)
                        , stat_sigma_(nbCluster)
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/

ParametersHandler<Clust::Gaussian_sk_>::ParametersHandler( ParametersHandler const& param)
                        : stat_mean_(param.stat_mean_)
                        , stat_sigma_(param.stat_sigma_)
{}
/* destructor */

ParametersHandler<Clust::Gaussian_sk_>::~ParametersHandler() {}
/* copy operator */

ParametersHandler<Clust::Gaussian_sk_>& ParametersHandler<Clust::Gaussian_sk_>::operator=( ParametersHandler const& other)
{
  stat_mean_ = other.stat_mean_;
  stat_sigma_ = other.stat_sigma_;
  return *this;
}
/* update statistics of the parameters. */

void ParametersHandler<Clust::Gaussian_sk_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].update(param.mean_[k]);
    stat_sigma_[k].update(param.sigma_[k]);
  }
}
/* Set the computed statistics */

void ParametersHandler<Clust::Gaussian_sk_>::setStatistics(Parameters& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    param.mean_[k] = stat_mean_[k].mean();
    stat_mean_[k].release();
    param.sigma_[k] = stat_sigma_[k].mean();
    stat_sigma_[k].release();
  }
}
/* Release the computed statistics */

void ParametersHandler<Clust::Gaussian_sk_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    stat_mean_[k].release();
    stat_sigma_[k].release();
  }
}
/* Initialize the statistics. */

void ParametersHandler<Clust::Gaussian_sk_>::resize(Range const& range)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].resize(range);
    stat_sigma_[k].release();
  }
}

/* default constructor. */

ParametersHandler<Clust::Gaussian_sj_>::ParametersHandler( int nbCluster)
                        : stat_mean_(nbCluster)
                        , stat_sigma_()
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/

ParametersHandler<Clust::Gaussian_sj_>::ParametersHandler( ParametersHandler const& param)
                        : stat_mean_(param.stat_mean_)
                        , stat_sigma_(param.stat_sigma_)
{}
/* destructor */

ParametersHandler<Clust::Gaussian_sj_>::~ParametersHandler() {}
/* copy operator */

ParametersHandler<Clust::Gaussian_sj_>& ParametersHandler<Clust::Gaussian_sj_>::operator=( ParametersHandler const& other)
{
  stat_mean_ = other.stat_mean_;
  stat_sigma_ = other.stat_sigma_;
  return *this;
}
/* update statistics of the parameters. */

void ParametersHandler<Clust::Gaussian_sj_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].update(param.mean_[k]);}
  stat_sigma_.update(param.sigma_);
}
/* Set the computed statistics */

void ParametersHandler<Clust::Gaussian_sj_>::setStatistics(Parameters& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    param.mean_[k] = stat_mean_[k].mean();
    stat_mean_[k].release();
  }
  param.sigma_ = stat_sigma_.mean();
  stat_sigma_.release();
}
/* Release the computed statistics */

void ParametersHandler<Clust::Gaussian_sj_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].release();}
  stat_sigma_.release();
}
/* Initialize the statistics. */

void ParametersHandler<Clust::Gaussian_sj_>::resize(Range const& range)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].resize(range);}
  stat_sigma_.resize(range);
}


/* default constructor. */

ParametersHandler<Clust::Gaussian_s_>::ParametersHandler( int nbCluster)
                        : stat_mean_(nbCluster)
                        , stat_sigma_()
{}
/* copy constructor.
 *  @param param the parameters to copy.
 **/

ParametersHandler<Clust::Gaussian_s_>::ParametersHandler( ParametersHandler const& param)
                        : stat_mean_(param.stat_mean_)
                        , stat_sigma_(param.stat_sigma_)
{}
/* destructor */

ParametersHandler<Clust::Gaussian_s_>::~ParametersHandler() {}
/* copy operator */

ParametersHandler<Clust::Gaussian_s_>& ParametersHandler<Clust::Gaussian_s_>::operator=( ParametersHandler const& other)
{
  stat_mean_ = other.stat_mean_;
  stat_sigma_ = other.stat_sigma_;
  return *this;
}
/* update statistics of the parameters. */

void ParametersHandler<Clust::Gaussian_s_>::updateStatistics(Parameters const& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].update(param.mean_[k]);}
  stat_sigma_.update(param.sigma_);
}
/* Set the computed statistics */

void ParametersHandler<Clust::Gaussian_s_>::setStatistics(Parameters& param)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  {
    param.mean_[k] = stat_mean_[k].mean();
    stat_mean_[k].release();
  }
  param.sigma_ = stat_sigma_.mean();
  stat_sigma_.release();
}
/* Release the computed statistics */

void ParametersHandler<Clust::Gaussian_s_>::releaseStatistics()
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].release();}
  stat_sigma_.release();
}
/* Initialize the statistics. */

void ParametersHandler<Clust::Gaussian_s_>::resize(Range const& range)
{
  for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
  { stat_mean_[k].resize(range);}
}

} // namespace STK

