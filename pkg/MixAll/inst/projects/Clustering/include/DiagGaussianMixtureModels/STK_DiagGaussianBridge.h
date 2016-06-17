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

/** @file STK_DiagGaussianBridge.h
 *  @brief In this file we define the bridge classes between the diagonal
 *  Gaussian mixtures and the composer.
 **/

#ifndef STK_DIAGGAUSSIANBRIDGE_H
#define STK_DIAGGAUSSIANBRIDGE_H

#include "STK_MixtureGaussian_s.h"
#include "STK_MixtureGaussian_sj.h"
#include "STK_MixtureGaussian_sjk.h"
#include "STK_MixtureGaussian_sk.h"

#include "../STK_IMixtureBridge.h"

#include <DManager/include/STK_DataBridge.h>
#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{

// forward declaration
template<int Id, class Data> class DiagGaussianBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixtureGaussian_sjk model
 **/
template<class Data_>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_sjk_, Data_> >
{
  typedef Data_ Data;
  /** Type of the mixture model */
  typedef MixtureGaussian_sjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gaussian_sjk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gaussian_sjk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixtureGaussian_sk model
 **/
template<class Data_>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_sk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGaussian_sk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gaussian_sk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gaussian_sk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixtureGaussian_sj model
 **/
template<class Data_>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_sj_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGaussian_sj<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gaussian_sj_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gaussian_sj_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixtureGaussian_s model
 **/
template<class Data_>
struct MixtureBridgeTraits< DiagGaussianBridge< Clust::Gaussian_s_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureGaussian_s<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Gaussian_s_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Gaussian_s_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Gaussian_
  };
};

} // namespace hidden

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGaussian_sjk model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_sjk_>
{
    typedef ModelParameters<Clust::Gaussian_sjk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_sigma_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_mean_(nbCluster)
                            , stat_sigma_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_mean_(param.stat_mean_)
                            , stat_sigma_(param.stat_sigma_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_mean_ = other.stat_mean_;
      stat_sigma_ = other.stat_sigma_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].update(param.mean_[k]);
        stat_sigma_[k].update(param.sigma_[k]);
      }
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      {
        param.mean_[k] = stat_mean_[k].mean();
        stat_mean_[k].release();
        param.sigma_[k] = stat_sigma_[k].mean();
        stat_sigma_[k].release();
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
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.mean_[k];
        params.row(kp+1) = param.sigma_[k];
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
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        param.mean_[k] = params.row(kp);
        param.sigma_[k] = params.row(kp+1);
      }
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      {
        stat_mean_[k].release();
        stat_sigma_[k].release();
      }
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].resize(range);
        stat_sigma_[k].resize(range);
      }
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGaussian_sk model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_sk_>
{
    typedef ModelParameters<Clust::Gaussian_sk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_sigma_;

    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_mean_(nbCluster)
                            , stat_sigma_(nbCluster)
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_mean_(param.stat_mean_)
                            , stat_sigma_(param.stat_sigma_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_mean_ = other.stat_mean_;
      stat_sigma_ = other.stat_sigma_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].update(param.mean_[k]);
        stat_sigma_[k].update(param.sigma_[k]);
      }
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      {
        param.mean_[k] = stat_mean_[k].mean();
        stat_mean_[k].release();
        param.sigma_[k] = stat_sigma_[k].mean();
        stat_sigma_[k].release();
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
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.mean_[k];
        params.row(kp+1) = param.sigma_[k];
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
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        param.mean_[k] = params.row(kp);
        param.sigma_[k] = params.row(kp+1).mean();
      }
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      {
        stat_mean_[k].release();
        stat_sigma_[k].release();
      }
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].resize(range);
        stat_sigma_[k].release();
      }
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGaussian_sj model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_sj_>
{
    typedef ModelParameters<Clust::Gaussian_sj_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
   Stat::Online<CPointX, Real> stat_sigma_;

    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_mean_(nbCluster)
                            , stat_sigma_()
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_mean_(param.stat_mean_)
                            , stat_sigma_(param.stat_sigma_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_mean_ = other.stat_mean_;
      stat_sigma_ = other.stat_sigma_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].update(param.mean_[k]);}
      stat_sigma_.update(param.sigma_);
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
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.mean_[k];
        params.row(kp+1) = param.sigma_;
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
      param.sigma_ =0.;
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        param.mean_[k] = params.row(kp);
        param.sigma_ += params.row(kp+1);
      }
      param.sigma_ /= param.mean_.size();
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      {
        param.mean_[k] = stat_mean_[k].mean();
        stat_mean_[k].release();
      }
      param.sigma_ = stat_sigma_.mean();
      stat_sigma_.release();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].release();}
      stat_sigma_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].resize(range);}
      stat_sigma_.resize(range);
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixtureGaussian_s model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_s_>
{
    typedef ModelParameters<Clust::Gaussian_s_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_sigma_;
    /** default constructor. */
    inline ParametersHandler( int nbCluster)
                            : stat_mean_(nbCluster)
                            , stat_sigma_()
    {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_mean_(param.stat_mean_)
                            , stat_sigma_(param.stat_sigma_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_mean_ = other.stat_mean_;
      stat_sigma_ = other.stat_sigma_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].update(param.mean_[k]);}
      stat_sigma_.update(param.sigma_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      {
        param.mean_[k] = stat_mean_[k].mean();
        stat_mean_[k].release();
      }
      param.sigma_ = stat_sigma_.mean();
      stat_sigma_.release();
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
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.mean_[k];
        params.row(kp+1) = param.sigma_;
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
      param.sigma_ =0.;
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        param.mean_[k] = params.row(kp);
        param.sigma_ += params.row(kp+1).mean();
      }
      param.sigma_ /= param.mean_.size();
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].release();}
      stat_sigma_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_mean_.begin(); k<stat_mean_.end(); ++k)
      { stat_mean_[k].resize(range);}
    }
};

/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 *  This class inherit from the interface IMixture and delegate almost
 *  all the treatments to the wrapped class.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 */
template<int Id, class Data>
class DiagGaussianBridge: public IMixtureBridge< DiagGaussianBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< DiagGaussianBridge<Id,Data> > Base;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits< DiagGaussianBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< DiagGaussianBridge<Id,Data> >::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits< DiagGaussianBridge<Id,Data> >::ParamHandler ParamHandler;
    // type of data
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Gaussian_
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
    DiagGaussianBridge( DataBridge<Data>* p_data, std::string const& idData, int nbCluster)
                      : Base( p_data, idData, nbCluster)
    { removeMissing(); initializeBridge();}
    /** copy constructor */
    DiagGaussianBridge( DiagGaussianBridge const& bridge): Base(bridge)
    { initializeBridge();}
    /** destructor */
    virtual ~DiagGaussianBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual DiagGaussianBridge* clone() const { return new DiagGaussianBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual DiagGaussianBridge* create() const
    {
      DiagGaussianBridge* p_bridge = new DiagGaussianBridge( mixture_, this->idData(), this->nbCluster());
      p_bridge->p_data_ = p_data_;
      p_bridge->initializeBridge();
      return p_bridge;
    }
    /** This function is used in order to get the current values of the
     *  parameters in an array.
     *  @param[out] params the array with the parameters of the mixture.
     */
    template<class Array>
    void getParameters(Array& params) const;
    /** Write the parameters on the output stream os */
    void writeParameters(ostream& os) const;

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
    DiagGaussianBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                      : Base(mixture, idData, nbCluster)
    {}
};

// implementation
template<int Id, class Data>
void DiagGaussianBridge<Id, Data>::removeMissing()
{
  Type value = Type();
  int j, old_j = Arithmetic<int>::NA();
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  {
    j = it->second; // get column
    if (j != old_j)
    {
      old_j = j;
      value =  p_data_->dataij().col(j).meanSafe();
    }
    p_data_->dataij()(it->first, j) = value;
  }
}

template<int Id, class Data>
template<class Array>
void DiagGaussianBridge<Id, Data>::getParameters(Array& params) const
{
  int nbClust = this->nbCluster();
  params.resize(2*nbClust, mixture_.p_data()->cols());
  for (int k= 0; k < nbClust; ++k)
  {
    for (int j= params.beginCols();  j< params.endCols(); ++j)
    {
      params(baseIdx+2*k  , j) = mixture_.mean(baseIdx + k, j);
      params(baseIdx+2*k+1, j) = mixture_.sigma(baseIdx + k, j);
    }
  }
}

/** Write the parameters on the output stream os */
template<int Id, class Data>
void DiagGaussianBridge<Id, Data>::writeParameters(ostream& os) const
{
  PointX m(mixture_.p_data()->cols());
  PointX s(mixture_.p_data()->cols());
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    // store sigma values in an array for a nice output
    for (int j= s.begin();  j < s.end(); ++j)
    { m[j] = mixture_.mean(k,j); s[j] = mixture_.sigma(k,j);}
    os << _T("---> Component ") << k << _T("\n");
    os << _T("mean = ") << m;
    os << _T("sigma = ")<< s;
  }
}


} // namespace STK

#endif /* STK_DIAGGAUSSIANBRIDGE_H */
