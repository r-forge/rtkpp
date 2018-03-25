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

/** @file STK_HDGaussianBridge.h
 *  @brief In this file we define the bridge classes between the diagonal
 *  Gaussian mixtures and the composer.
 **/

#ifndef STK_HDGAUSSIANBRIDGE_H
#define STK_HDGAUSSIANBRIDGE_H

//#include "STK_DiagGaussian_s.h"
//#include "STK_DiagGaussian_sj.h"
//#include "STK_DiagGaussian_sjk.h"
//#include "STK_DiagGaussian_sk.h"

#include "../STK_IMixtureBridge.h"

#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{

// forward declaration
template<int Id, class Data> class HDGaussianBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the HDGaussian_ajk_bk_qk_d_ model
 **/
template<class Data_>
struct MixtureBridgeTraits< HDGaussianBridge< Clust::HDGaussian_ajk_bk_qk_d_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef HDGaussian_ajk_bk_qk_d<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::HDGaussian_ajk_bk_qk_d_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::HDGaussian_ajk_bk_qk_d_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::HDGaussian_
  };
};

} // namespace hidden

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for DiagGaussian_sjk model
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
    ParametersHandler( int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ParametersHandler( ParametersHandler const& param);
    /** destructor */
    ~ParametersHandler();
    /** copy operator */
    ParametersHandler& operator=( ParametersHandler const& other);
    /** update statistics of the parameters. */
    void updateStatistics(Parameters const& param);
    /** Set the computed statistics */
    void setStatistics(Parameters& param);
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the mean and
     *  sd parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        params.row(kp) = param.mean_[k];
        params.row(kp+1) = param.sigma_[k];
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the mean and
     *  sd parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      for(int k=param.mean_.begin(), kp= params.beginRows(); k<param.mean_.end(); ++k, kp+=2)
      {
        param.mean_[k] = params.row(kp);
        param.sigma_[k] = params.row(kp+1);
      }
    }
    /** Release the computed statistics */
    void releaseStatistics();
    /** Initialize the statistics. */
    void resize(Range const& range);
};


/** @ingroup Clustering
 *  @brief template implementation of the IMixture interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 *  This class inherit from the interface IMixture and delegate almost
 *  all the treatments to the wrapped class. The bridge handles the missing
 *  values and the averaging of the parameters and imputed/simulated mising
 *  values during the estimation process.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class. This identifier should be find in the
 * Clust::Mixture enum.
 * @tparam Data container of the data used by the STK::DataBridge class
 */
template<int Id, class Data>
class HDGaussianBridge: public IMixtureBridge< HDGaussianBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< HDGaussianBridge<Id,Data> > Base;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits< HDGaussianBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< HDGaussianBridge<Id,Data> >::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits< HDGaussianBridge<Id,Data> >::ParamHandler ParamHandler;
    // type of data
    typedef typename hidden::MixtureBridgeTraits< HDGaussianBridge<Id,Data> >::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::HDGaussian_
    };
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::paramHandler_;
    using Base::p_dataij_;
    using Base::p_tik;
    using Base::v_missing_;
    using Base::removeMissing;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_dataij pointer on the data set used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    HDGaussianBridge( Data* p_dataij, std::string const& idData, int nbCluster)
                    : Base( p_dataij, idData, nbCluster)
    {
      removeMissing(); // remove missing from data only once at creation
      mixture_.setData(p_data_->dataij());
      paramHandler_.resize(p_data_->cols());
    }
    /** copy constructor */
    HDGaussianBridge( HDGaussianBridge const& bridge): Base(bridge)
    {
      mixture_.setData(p_data_->dataij());
      paramHandler_.resize(p_data_->cols());
    }
    /** destructor */
    virtual ~HDGaussianBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual HDGaussianBridge* clone() const { return new HDGaussianBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual HDGaussianBridge* create() const
    {
      HDGaussianBridge* p_bridge = new HDGaussianBridge( mixture_, this->idData(), this->nbCluster());
      p_bridge->p_data_ = p_data_;
      p_bridge->p_dataij_ = p_dataij_;
      p_bridge->mixture_.setData(*p_dataij_);
      p_bridge->paramHandler_.resize(p_data_->cols());
      p_bridge->v_missing_ = v_missing_;
      return p_bridge;
    }
    /** @return a safe value for the jth variable
     *  @param j index of the column with the safe value needed */
    Type safeValue( int j) const
    { return p_data_->dataij().col(j).meanSafe();}

  private:
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    HDGaussianBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                   : Base(mixture, idData, nbCluster)
    {}
};

} // namespace STK

#endif /* STK_DIAGGAUSSIANBRIDGE_H */
