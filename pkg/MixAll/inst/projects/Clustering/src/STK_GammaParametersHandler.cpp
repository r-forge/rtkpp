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

/** @file STK_GammaParametersHandler.cpp
 *  @brief In this file we implement the parameters handlers classes for categorical
 *  mixture models
 **/

#include "../include/GammaModels/STK_GammaParametersHandler.h"

namespace STK
{
    /* default constructor. */
    ParametersHandler<Clust::Gamma_ajk_bjk_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ajk_bjk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ajk_bjk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ajk_bjk_>& ParametersHandler<Clust::Gamma_ajk_bjk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ajk_bjk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_bjk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_bjk_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ajk_bjk_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].resize(range);
        stat_scale_[k].resize(range);
      }
    }

    /* default constructor. */
    ParametersHandler<Clust::Gamma_ajk_bk_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ajk_bk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ajk_bk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ajk_bk_>& ParametersHandler<Clust::Gamma_ajk_bk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ajk_bk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_bk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_bk_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ajk_bk_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].resize(range);
        stat_scale_[k].release();
      }
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_ajk_bj model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_ajk_bj_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ajk_bj_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ajk_bj_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ajk_bj_>& ParametersHandler<Clust::Gamma_ajk_bj_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ajk_bj_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
      }
      stat_scale_.update(param.scale_);
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_bj_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_bj_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ajk_bj_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].resize(range);
      }
      stat_scale_.resize(range);
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_ajk_b model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_ajk_b_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ajk_b_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ajk_b_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ajk_b_>& ParametersHandler<Clust::Gamma_ajk_b_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ajk_b_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);}
      stat_scale_.update(param.scale_);
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_b_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_ajk_b_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ajk_b_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].resize(range);
      }
      stat_scale_.release();
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_ak_bjk model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_ak_bjk_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ak_bjk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ak_bjk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ak_bjk_>& ParametersHandler<Clust::Gamma_ak_bjk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ak_bjk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_bjk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_bjk_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ak_bjk_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();
        stat_scale_[k].resize(range);
      }
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_ak_bk model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_ak_bk_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ak_bk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ak_bk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ak_bk_>& ParametersHandler<Clust::Gamma_ak_bk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ak_bk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_bk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_bk_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ak_bk_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_ak_bj model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_ak_bj_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ak_bj_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ak_bj_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ak_bj_>& ParametersHandler<Clust::Gamma_ak_bj_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ak_bj_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
      }
      stat_scale_.update(param.scale_);
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_bj_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_bj_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ak_bj_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
      }
      stat_scale_.resize(range);
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_ak_b model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_ak_b_>::ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_ak_b_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_ak_b_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_ak_b_>& ParametersHandler<Clust::Gamma_ak_b_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_ak_b_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
      }
      stat_scale_.update(param.scale_);
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_b_>::setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_ak_b_>::releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_ak_b_>::resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
      }
      stat_scale_.release();
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_aj_bjk model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_aj_bjk_>::ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_aj_bjk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_aj_bjk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_aj_bjk_>& ParametersHandler<Clust::Gamma_aj_bjk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_aj_bjk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
      stat_shape_.update(param.shape_);

    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_aj_bjk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_aj_bjk_>::releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_aj_bjk_>::resize(Range const& range)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].resize(range);}
      stat_shape_.resize(range);
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_aj_bk model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_aj_bk_>::ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_aj_bk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_aj_bk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_aj_bk_>& ParametersHandler<Clust::Gamma_aj_bk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_aj_bk_>::updateStatistics(Parameters const& param)
    {
      stat_shape_.update(param.shape_);
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_aj_bk_>::setStatistics(Parameters& param)
    {
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_aj_bk_>::releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_aj_bk_>::resize(Range const& range)
    {
      stat_shape_.resize(range);
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_aj_bjk model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_a_bjk_>::ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_a_bjk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_a_bjk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_a_bjk_>& ParametersHandler<Clust::Gamma_a_bjk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_a_bjk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
      stat_shape_.update(param.shape_);
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_a_bjk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
    }
    /* Release the computed statistics */
    void ParametersHandler<Clust::Gamma_a_bjk_>::releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_a_bjk_>::resize(Range const& range)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].resize(range);}
      stat_shape_.release();
    }

/* @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_a_bk model
 **/
    /* default constructor. */
    ParametersHandler<Clust::Gamma_a_bk_>::ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /* copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler<Clust::Gamma_a_bk_>::ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /* destructor */
    ParametersHandler<Clust::Gamma_a_bk_>::~ParametersHandler() {}
    /* copy operator */
    ParametersHandler<Clust::Gamma_a_bk_>& ParametersHandler<Clust::Gamma_a_bk_>::operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /* update statistics of the parameters. */
    void ParametersHandler<Clust::Gamma_a_bk_>::updateStatistics(Parameters const& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
      stat_shape_.update(param.shape_);
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_a_bk_>::setStatistics(Parameters& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
    }
    /* Set the computed statistics */
    void ParametersHandler<Clust::Gamma_a_bk_>::releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /* Initialize the statistics. */
    void ParametersHandler<Clust::Gamma_a_bk_>::resize(Range const& range)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }


} // namespace STK
