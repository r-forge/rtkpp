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
 * created on: Dec 4, 2013
 * Authors: Serge Iovleff
 **/

/** @file STK_MixtureGammaBase.h
 *  @brief In this file we implement the base class for the gamma models
 **/

#ifndef STK_MIXTUREGAMMABASE_H
#define STK_MIXTUREGAMMABASE_H

#include "../STK_IMixtureDensity.h"

#include <STatistiK/include/STK_Law_Categorical.h>
#include <STatistiK/include/STK_Law_Gamma.h>
#include <STatistiK/include/STK_Stat_Functors.h>

#include <Analysis/include/STK_Algo_FindZero.h>
#include <Analysis/include/STK_Funct_raw.h>
#include <Analysis/include/STK_Funct_gamma.h>

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  Functor computing the derivative of the lnLikelihood of a gamma_ajk_bjk model */
class invPsiMLog : public IFunction<invPsiMLog >
{
  public:
    inline invPsiMLog( Real const& y): y_(y)  {}
    /** @return the value of the function at a
     * @param a a positive real value
     **/
    inline Real fImpl(Real const& a) const
    { return (y_ + std::log(a) - Funct::psi_raw(a));}
    /** @return the minimal value of the function at x */
    inline Real xminImpl() const { return 0;}

  private:
    Real y_;
};

/** @ingroup hidden
 *  Functor computing the difference between the psi function  and a fixed value
 **/
class invPsi : public IFunction<invPsi >
{
  public:
    /** initialize y_ */
    inline invPsi( Real const& y): y_(y) {}
    /** @return the value of the function at a
     *  @param x a positive real value
     **/
    inline Real fImpl(Real const& x) const { return (y_ - Funct::psi_raw(x));}
    /** @return the minimal value of the function at x */
    inline Real xminImpl() const { return 0;}
  private:
    Real y_;
};

} // namespace hidden

/** base class of the Gamma models Parameter Handler */
struct ParametersGammaBase
{
  /** mean */
  Array1D<CPointX> mean_;
  /** log-means */
  Array1D<CPointX> meanLog_;
  /** variance_ */
  Array1D<CPointX> variance_;
  /** copy operator */
  inline ParametersGammaBase& operator=( ParametersGammaBase const& other)
  { mean_ = other.mean_;
    meanLog_ = other.meanLog_;
    variance_ = other.variance_;
    return *this;
  }

  /** default constructor */
  ParametersGammaBase( int nbCluster)
                     : mean_(nbCluster)
                     , meanLog_(nbCluster)
                     , variance_(nbCluster) {}
  /** copy constructor */
  ParametersGammaBase( ParametersGammaBase const& model)
                     : mean_(model.mean_)
                     , meanLog_(model.meanLog_)
                     , variance_(model.variance_) {}
  /** destructor */
  inline ~ParametersGammaBase() {}
  /** @param range the range of the variables */
  void resize(Range const& range)
  {
    for (int k=mean_.begin(); k < mean_.end(); ++k)
    {
      mean_[k].resize(range) = 1.;
      meanLog_[k].resize(range) = 0.;
      variance_[k].resize(range) = 1.;
    }
  }
};

/** @ingroup Clustering
 *  Base class for the gamma models
 **/
template<class Derived>
class MixtureGammaBase : public IMixtureDensity<Derived >
{
  public:
    typedef IMixtureDensity<Derived > Base;
    using Base::param_;
    using Base::p_data;

  protected:
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline MixtureGammaBase( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixtureGammaBase( MixtureGammaBase const& model): Base(model) {}
    /** destructor */
    inline ~MixtureGammaBase() {}

  public:
    /** @return the shape of the kth cluster and jth variable */
    inline Real shape(int k, int j) const { return param_.shape(k,j);}
    /** @return the scale of the kth cluster and jth variable */
    inline Real scale(int k, int j) const { return param_.scale(k,j);}
    /** Initialize the parameters of the model. */
    void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param pk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    Real impute(int i, int j, Weights const& pk) const;
    /** @return a simulated value for the jth variable of the ith sample
     *  in the kth cluster.
     *  @param i,j,k indexes of the data to simulate
     **/
    inline Real rand(int i, int j, int k) const
    { return Law::Gamma::rand(shape(k,j), scale(k,j));}

  protected:
    /** compute the Q(theta) value. */
    Real qValue(CArrayXX const* const& p_tik, CPointX const* p_nk) const;
    /** compute the weighted moments of a gamma mixture. */
    bool moments(CArrayXX const* const& p_tik);
    /** get the weighted mean of the jth variable of the kth cluster. */
    inline Real meanjk( int j, int k) { return param_.mean_[k][j];}
    /** get the weighted variance of the jth variable of the kth cluster. */
    inline Real variancejk( int j, int k) { return param_.variance_[k][j];}
    /** get the mean of the weighted means of the kth cluster. */
    inline Real meank( int k) { return param_.mean_[k].mean();}
    /** get the mean of the weighted variances of the kth cluster. */
    inline Real variancek( int k) { return param_.variance_[k].mean();}
};

template<class Derived>
template<class Weights>
Real MixtureGammaBase<Derived>::impute(int i, int j, Weights const& pk) const
{
  Real sum = 0.;
  for (int k= pk.begin(); k < pk.end(); ++k)
  { sum += pk[k] * shape(k,j) * scale(k,j);}
  return sum;
}

/* compute safely the weighted moments of a gamma law. */
template<class Derived>
bool MixtureGammaBase<Derived>::moments(CArrayXX const* const& p_tik)
{
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    CVectorX tikColk(p_tik->col(k), true); // create a reference
    for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
    {
      // mean
      Real mean =  p_data()->col(j).wmean(tikColk);
      if ( (mean<=0) || isNA(mean) ) { return false;}
      param_.mean_[k][j] = mean;
      // mean log
      Real meanLog =  p_data()->col(j).log().wmean(tikColk);
      if (isNA(meanLog)) { return false;}
      param_.meanLog_[k][j] = meanLog;
      // variance
      Real variance =  p_data()->col(j).wvariance(mean, tikColk);
      if ((variance<=0)||isNA(variance)){ return false;}
      param_.variance_[k][j] = variance;
    }
  }
  return true;
}

/*get the parameters of the model*/
template<class Derived>
Real MixtureGammaBase<Derived>::qValue(CArrayXX const* const& p_tik, CPointX const* p_nk) const
{
  Real value = 0.;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    Real sumjk=0.0;
    for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
    {
      Real a = shape(k,j), b = scale(k,j);
      sumjk += a * (param_.meanLog_[k][j]-std::log(b))
             - param_.mean_[k][j]/b - STK::Funct::lgamma(a);
    }
    value += p_nk->elt(k) * sumjk;
  }
  return value;
}

} // namespace STK

#endif /* STK_MIXTUREGAMMABASE_H */
