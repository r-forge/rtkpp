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

/** @file STK_GammaBridge.h
 *  @brief In this file we define the bridge classes between the diagonal
 *  Gaussian mixtures and the composer.
 **/

#ifndef STK_GAMMABRIDGE_H
#define STK_GAMMABRIDGE_H

#include "STK_MixtureGamma_ajk_bjk.h"
#include "STK_MixtureGamma_ajk_bk.h"
#include "STK_MixtureGamma_ajk_bj.h"
#include "STK_MixtureGamma_ajk_b.h"
#include "STK_MixtureGamma_ak_bjk.h"
#include "STK_MixtureGamma_ak_bk.h"
#include "STK_MixtureGamma_ak_bj.h"
#include "STK_MixtureGamma_ak_b.h"
#include "STK_MixtureGamma_aj_bjk.h"
#include "STK_MixtureGamma_aj_bk.h"
#include "STK_MixtureGamma_a_bjk.h"
#include "STK_MixtureGamma_a_bk.h"

#include "../STK_IMixtureBridge.h"

#include <DManager/include/STK_DataBridge.h>
#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{

// forward declaration
template<int Id, class Data> class GammaBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ajk_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ajk_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_bjk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ajk_bjk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ajk_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ajk_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_bk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ajk_bk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ajk_bj_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_bj_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ajk_bj<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_bj_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ajk_bj_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ajk_b_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ajk_b_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ajk_b<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ajk_b_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ajk_b_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ak_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ak_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_bjk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ak_bjk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ak_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ak_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_bk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ak_bk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ak_bj_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_bj_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ak_bj<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_bj_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ak_bj_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_ak_b_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_ak_b_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_ak_b<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_ak_b_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_ak_b_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_aj_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_aj_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_aj_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_aj_bjk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_aj_bjk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_aj_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_aj_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_aj_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_aj_bk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_aj_bk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_a_bjk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_a_bjk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_a_bjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_a_bjk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_a_bjk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};
/** @ingroup hidden
 *  Partial  specialization of the MixtureBridgeTraits for the MixtureGamma_a_bk_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< GammaBridge< Clust::Gamma_a_bk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGamma_a_bk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gamma_a_bk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gamma_a_bk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gamma_
  };
};

} // namespace hidden

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ajk_bjk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_bjk_>
{
    typedef ModelParameters<Clust::Gamma_ajk_bjk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp);
        param.scale_[k] = params.row(kp+1);
      }
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].resize(range);
        stat_scale_[k].resize(range);
      }
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ajk_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_bk_>
{
    typedef ModelParameters<Clust::Gamma_ajk_bk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp);
        param.scale_[k] = params.row(kp+1).mean();
      }
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].resize(range);
        stat_scale_[k].release();
      }
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ajk_bj model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_bj_>
{
    typedef ModelParameters<Clust::Gamma_ajk_bj_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<CPointX, Real> stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
      }
      stat_scale_.update(param.scale_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_;
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.scale_ = 0.;
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp);
        param.scale_ += params.row(kp+1);
      }
      param.scale_ /= param.shape_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].resize(range);
      }
      stat_scale_.resize(range);
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ajk_b model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_b_>
{
    typedef ModelParameters<Clust::Gamma_ajk_b_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);}
      stat_scale_.update(param.scale_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_;
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.scale_ = 0.;
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp);
        param.scale_ = params.row(kp+1).mean();
      }
      param.scale_ /= param.shape_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].resize(range);
      }
      stat_scale_.release();
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ak_bjk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_bjk_>
{
    typedef ModelParameters<Clust::Gamma_ak_bjk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp).mean();
        param.scale_[k] = params.row(kp+1);
      }
    }
    /** Set the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();
        stat_scale_[k].resize(range);
      }
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ak_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_bk_>
{
    typedef ModelParameters<Clust::Gamma_ak_bk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
        stat_scale_[k].update(param.scale_[k]);
      }
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp).mean();
        param.scale_[k] = params.row(kp+1).mean();
      }
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();
        stat_scale_[k].release();
      }
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ak_bj model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_bj_>
{
    typedef ModelParameters<Clust::Gamma_ak_bj_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<CPointX, Real> stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
      }
      stat_scale_.update(param.scale_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_;
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.scale_ = 0.;
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp).mean();
        param.scale_ += params.row(kp+1);
      }
      param.scale_ /= param.shape_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
      }
      stat_scale_.resize(range);
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_ak_b model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_b_>
{
    typedef ModelParameters<Clust::Gamma_ak_b_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_(nbCluster)
                            , stat_scale_()
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].update(param.shape_[k]);
      }
      stat_scale_.update(param.scale_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        param.shape_[k] = stat_shape_[k].mean();
        stat_shape_[k].release();
      }
      param.scale_ = stat_scale_.mean();
      stat_scale_.release();
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_[k];
        params.row(kp+1) = param.scale_;
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.scale_ = 0.;
      for(int k=param.shape_.begin(), kp= params.beginRows(); k<param.shape_.end(); ++k, kp+=2)
      {
        param.shape_[k] = params.row(kp).mean();
        param.scale_ += params.row(kp+1).mean();
      }
      param.scale_ /= param.shape_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      { stat_shape_[k].release();}
      stat_scale_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_shape_.begin(); k<stat_shape_.end(); ++k)
      {
        stat_shape_[k].release();
      }
      stat_scale_.release();
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_aj_bjk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_aj_bjk_>
{
    typedef ModelParameters<Clust::Gamma_aj_bjk_> Parameters;
    /** shape statistics */
    Stat::Online<CPointX, Real> stat_shape_;
    /** Array of the scale statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
      stat_shape_.update(param.shape_);

    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_;
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.shape_ = 0.;
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        param.shape_ += params.row(kp);
        param.scale_[k] = params.row(kp+1);
      }
      param.shape_ /= param.scale_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].resize(range);}
      stat_shape_.resize(range);
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_aj_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_aj_bk_>
{
    typedef ModelParameters<Clust::Gamma_aj_bk_> Parameters;
    /** shape statistics */
    Stat::Online<CPointX, Real> stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      stat_shape_.update(param.shape_);
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_;
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.shape_ = 0.;
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        param.shape_ += params.row(kp);
        param.scale_[k] = params.row(kp+1).mean();
      }
      param.shape_ /= param.scale_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      stat_shape_.resize(range);
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_aj_bjk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_a_bjk_>
{
    typedef ModelParameters<Clust::Gamma_a_bjk_> Parameters;
    /** shape statistics */
    Stat::Online<Real, Real> stat_shape_;
    /** Array of the scale statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
      stat_shape_.update(param.shape_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_;
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.shape_ = 0.;
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        param.shape_ += params.row(kp).mean();
        param.scale_[k] = params.row(kp+1);
      }
      param.shape_ /= param.scale_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].resize(range);}
      stat_shape_.release();
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGamma_a_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_a_bk_>
{
    typedef ModelParameters<Clust::Gamma_a_bk_> Parameters;
    /** shape statistics */
    Stat::Online<Real, Real> stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_shape_()
                            , stat_scale_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_shape_(param.stat_shape_)
                            , stat_scale_(param.stat_scale_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_shape_ = other.stat_shape_;
      stat_scale_ = other.stat_scale_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].update(param.scale_[k]);}
      stat_shape_.update(param.shape_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      {
        param.scale_[k] = stat_scale_[k].mean();
        stat_scale_[k].release();
      }
      param.shape_ = stat_shape_.mean();
      stat_shape_.release();
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.shape_;
        params.row(kp+1) = param.scale_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      param.shape_ = 0.;
      for(int k=param.scale_.begin(), kp= params.beginRows(); k<param.scale_.end(); ++k, kp+=2)
      {
        param.shape_ += params.row(kp).mean();
        param.scale_[k] = params.row(kp+1).mean();
      }
      param.shape_ /= param.scale_.size();
    }
    /** Set the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_scale_.begin(); k<stat_scale_.end(); ++k)
      { stat_scale_[k].release();}
      stat_shape_.release();
    }
};

/** @ingroup Clustering
 *  @brief Templated implementation of the IMixtureBridge interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 *  This class inherit from the interface IMixture and delegate almost
 *  all the treatments to the wrapped class.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 */
template<int Id, class Data>
class GammaBridge: public IMixtureBridge< GammaBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< GammaBridge<Id,Data> > Base;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits< GammaBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< GammaBridge<Id,Data> >::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits< GammaBridge<Id,Data> >::ParamHandler ParamHandler;
    // type of data
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Gamma_
    };
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::paramHandler_;
    using Base::p_data_;
    using Base::p_tik;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_data pointer on the DataBridge that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    GammaBridge( DataBridge<Data>* p_data, std::string const& idData, int nbCluster)
               : Base( p_data, idData, nbCluster)
    { removeMissing(); initializeBridge();}
    /** copy constructor */
    GammaBridge( GammaBridge const& bridge): Base(bridge)
    { initializeBridge();}
    /** destructor */
    virtual ~GammaBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual GammaBridge* clone() const { return new GammaBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual GammaBridge* create() const
    {
      GammaBridge* p_bridge = new GammaBridge( mixture_, this->idData(), this->nbCluster());
      p_bridge->p_data_ = p_data_;
      p_bridge->initializeBridge();
      return p_bridge;
    }
    /** This function is used in order to get the current values of the parameters.
     *  @param params the array with the parameters of the mixture.
     */
    template<class Array>
    void getParameters(Array& params) const;
    /** This function can be used to write summary of parameters to the output stream.
     *  @param os Stream where you want to write the summary of parameters.
     */
    virtual void writeParameters(ostream& os) const;

  private:
    /** This function will be used for the imputation of the missing data
     *  at the initialization.
     **/
    void removeMissing();
    /** This function will be used in order to initialize the mixture model
     *  using informations stored by the DataBridge. For example the missing
     *  values in the case of a DataBridge instance.
     **/
    void initializeBridge()
    {
      mixture_.setData(p_data_->dataij());
      paramHandler_.resize(p_data_->cols());
    }
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    GammaBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
               : Base(mixture, idData, nbCluster)
    {}
};

// implementation
template<int Id, class Data>
void GammaBridge<Id, Data>::removeMissing()
{
  Type value = Type();
  int j, old_j = Arithmetic<int>::NA();
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  {
    j = it->second; // get column
    if (j != old_j)
    {
      old_j = j;
      value =  p_data_->dataij().col(j).safe(1).mean();
    }
    p_data_->dataij()(it->first, j) = value;
  }
}

template<int Id, class Data>
template<class Array>
void GammaBridge<Id, Data>::getParameters(Array& params) const
{
  int nbClust = this->nbCluster();
  params.resize(2*nbClust, mixture_.p_data()->cols());
  for (int k= 0; k < nbClust; ++k)
  {
    for (int j=  mixture_.p_data()->beginCols();  j< mixture_.p_data()->endCols(); ++j)
    {
      params(baseIdx+2*k  , j) = mixture_.shape(baseIdx+k,j);
      params(baseIdx+2*k+1, j) = mixture_.scale(baseIdx+k,j);
    }
  }
}

template<int Id, class Data>
void GammaBridge<Id, Data>::writeParameters(ostream& os) const
{
    Array2DPoint<Real> a(mixture_.p_data()->cols()), b(mixture_.p_data()->cols());
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    {
      // store shape and scale values in an array for a nice output
      for (int j=mixture_.p_data()->beginCols();  j < mixture_.p_data()->endCols(); ++j)
      {
        a[j] = mixture_.shape(k,j);
        b[j] = mixture_.scale(k,j);
      }
      os << _T("---> Component ") << k << _T("\n");
      os << _T("shape = ") << a;
      os << _T("scale = ") << b;
    }
}


} // namespace STK

#endif /* STK_GAMMABRIDGE_H */
