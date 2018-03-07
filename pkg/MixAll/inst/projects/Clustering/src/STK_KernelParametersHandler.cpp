/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * Project:  stkpp::
 * created on: Sep 8, 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_KernelParametersHandler.cpp
 *  @brief In this file we implement the parameters handlers classes for kernel
 *  mixture models.
 **/


#include "../include/KernelMixtureModels/STK_KernelParametersHandler.h"

namespace STK
{
/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for KernelGaussian_s model
 **/
    /* default constructor. All lambdas are initialized to 1. */
    ParametersHandler<Clust::KernelGaussian_s_>::ParametersHandler( int nbCluster)
                            : stat_sigma2_(), stat_dim_(nbCluster)
    {}
    /* copy constructor.
     * @param param the parameters to copy.
     **/
    ParametersHandler<Clust::KernelGaussian_s_>::ParametersHandler( ParametersHandler const& param)
                            : stat_sigma2_(param.stat_sigma2_), stat_dim_(param.stat_dim_)
    {}
    /* destructor */
    ParametersHandler<Clust::KernelGaussian_s_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::KernelGaussian_s_>& ParametersHandler<Clust::KernelGaussian_s_>::operator=( ParametersHandler const& other)
    {
      stat_sigma2_ = other.stat_sigma2_;
      stat_dim_ = other.stat_dim_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::KernelGaussian_s_>::updateStatistics(Parameters const& param)
    {
      stat_sigma2_.update(param.sigma2_);
      for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
      { stat_dim_[k].update(param.dim_[k]);}
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::KernelGaussian_s_>::setStatistics(Parameters& param)
    {
      param.sigma2_ = stat_sigma2_.mean();
      stat_sigma2_.release();
      for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
      {
        param.dim_[k] = stat_dim_[k].mean();
        stat_dim_[k].release();
      }
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::KernelGaussian_s_>::releaseStatistics()
    {
      stat_sigma2_.release();
      for(int k=stat_dim_.begin(); k<stat_dim_.end(); ++k)
      { stat_dim_[k].release();}
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for KernelGaussian_sk model
 **/
    /* default constructor. All lambdas are initialized to 1. */
    ParametersHandler<Clust::KernelGaussian_sk_>::ParametersHandler( int nbCluster)
                            : stat_sigma2_(nbCluster), stat_dim_(nbCluster)
    {}
    /* copy constructor.
     * @param param the parameters to copy.
     **/
    ParametersHandler<Clust::KernelGaussian_sk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_sigma2_(param.stat_sigma2_), stat_dim_(param.stat_dim_)
    {}
    /* destructor */
    ParametersHandler<Clust::KernelGaussian_sk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::KernelGaussian_sk_>& ParametersHandler<Clust::KernelGaussian_sk_>::operator=( ParametersHandler const& other)
    {
      stat_sigma2_ = other.stat_sigma2_;
      stat_dim_ = other.stat_dim_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::KernelGaussian_sk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_sigma2_.begin(); k<stat_sigma2_.end(); ++k)
      { stat_sigma2_[k].update(param.sigma2_[k]);
        stat_dim_[k].update(param.dim_[k]);
      }
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::KernelGaussian_sk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_sigma2_.begin(); k<stat_sigma2_.end(); ++k)
      {
        param.sigma2_[k] = stat_sigma2_[k].mean();
        stat_sigma2_[k].release();
        param.dim_[k] = stat_dim_[k].mean();
        stat_dim_[k].release();
      }
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::KernelGaussian_sk_>::releaseStatistics()
    {
      for(int k=stat_sigma2_.begin(); k<stat_sigma2_.end(); ++k)
      {
        stat_sigma2_[k].release();
        stat_dim_[k].release();
      }
    }

} // namespace STK
