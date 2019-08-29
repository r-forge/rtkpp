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

/** @file STK_HDGaussianBase.h
 *  @brief In this file we implement the base class for the HD Gaussian models
 **/

#ifndef STK_HDGAUSSIANBASE_H
#define STK_HDGAUSSIANBASE_H

#include "../STK_IMixtureDensity.h"
#include <Arrays/include/STK_CArray.h>
#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <Arrays/include/STK_Const_Arrays.h>
#include <Arrays/include/STK_Display.h>
#include <STatistiK/include/STK_Law_Normal.h>
#include <STatistiK/include/STK_Law_UniformDiscrete.h>
#include <STatistiK/include/STK_Stat_Functors.h>

namespace STK
{

/** @ingroup Clustering
 *  Base class for the diagonal Gaussian models
 **/
template<class Derived>
class HDGaussianBase: public IMixtureDensity<Derived >
{
  public:
    typedef IMixtureDensity<Derived > Base;
    using Base::nbCluster;
    using Base::param_;
    using Base::p_data;

  protected:
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline HDGaussianBase( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline HDGaussianBase( HDGaussianBase const& model): Base(model) {}
    /** destructor */
    inline ~HDGaussianBase() {}

  public:
    /** @return mean of the kth cluster and jth variable */
    inline Real const& mean(int k, int j) const { return param_.mean(k,j);}
    /** @return dimension of low sub-subspace of the kth cluster */
    inline int const& d(int k) const { return param_.d(k);}
    /** @return noise variance of the kth cluster */
    inline Real const& b(int k) const { return param_.b(k);}
    /** @return rotation matrix of the kth cluster */
    inline Real const& q(int k) const { return param_.q(k);}
    /** @return rotation matrix of the kth cluster */
    inline Real const& a(int j, int k) const { return param_.a(k, j);}

    /** Initialize the parameters of the model. */
    inline void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param pk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    Real impute(int i, int j, Weights const& pk) const;
    /** @return a simulated value for the jth variable of the ith sample
     * in the kth cluster
     * @param i,j,k indexes of the data to simulate */
    inline Real rand(int i, int j, int k) const
    { return Law::Normal::rand(mean(k, j), sigma(k,j));}
    /** This function is used in order to get the current values of the means
     *  and standard deviations.
     *  @param[out] params the array with the parameters of the mixture.
     */
    template<class Array>
    void getParameters(Array& params) const;
    /** This function can be used to write summary of parameters to the output stream.
     *  @param p_tik a constant pointer on the posterior probabilities
     *  @param os Stream where you want to write the summary of parameters.
     */
    void writeParameters(CArrayXX const* p_tik, ostream& os) const;

  protected:
    /** sample randomly the mean of each component by sampling randomly a row
     *  of the data set.
     **/
    void randomMean( CArrayXX const*  p_tik);
    /** compute the weighted mean of a Gaussian mixture. */
    bool updateMean( CArrayXX const*  p_tik);
};

template<class Derived>
template<class Weights>
Real HDGaussianBase<Derived>::impute(int i, int j, Weights const& pk) const
{
  Real sum = 0.;
  for (int k= pk.begin(); k < pk.end(); ++k)
  { sum += pk[k] * mean(k,j);}
  return sum;
}

template<class Derived>
void HDGaussianBase<Derived>::randomMean( CArrayXX const*  p_tik)
{
  // indexes array
  VectorXi indexes(this->nbSample());
  for(int i=p_data()->beginRows(); i< p_data()->endRows(); ++i) { indexes[i] = i;}
  Range rind = this->nbSample();
  // sample without repetition
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    // random number in [0, end-k[
    int i = Law::UniformDiscrete::rand(rind.begin(), rind.end()-1);
    // get ith individuals
    param_.mean_[k] = p_data()->row(indexes[i]);
    // exchange it with nth
    indexes.swap(i, rind.lastIdx());
    // decrease
    rind.decLast(1);
  }
}

template<class Derived>
bool HDGaussianBase<Derived>::updateMean( CArrayXX const*  p_tik)
{
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
    { param_.mean_[k][j] = p_data()->col(j).wmean(p_tik->col(k));}
  }
  return true;
}

/* This function is used in order to get the current values of the means
 *  and standard deviations.
 *  @param[out] params the array with the parameters of the mixture.
 */
template<class Derived>
template<class Array>
void HDGaussianBase<Derived>::getParameters(Array& params) const
{
  int nbClust = nbCluster();
  params.resize(2*nbClust, p_data()->cols());
  for (int k= 0; k < nbClust; ++k)
  {
    for (int j= params.beginCols();  j< params.endCols(); ++j)
    {
      params(baseIdx+2*k  , j) = mean(baseIdx + k, j);
      params(baseIdx+2*k+1, j) = sigma(baseIdx + k, j);
    }
  }
}


/* This function can be used to write summary of parameters to the output stream.
 *  @param p_tik a constant pointer on the posterior probabilities
 *  @param os Stream where you want to write the summary of parameters.
 */
template<class Derived>
void HDGaussianBase<Derived>::writeParameters(CArrayXX const* p_tik, ostream& os) const
{
  CPointX m(p_data()->cols()), s(p_data()->cols());
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    // store sigma values in an array for a nice output
    for (int j= s.begin();  j < s.end(); ++j)
    { m[j] = mean(k,j); s[j] = sigma(k,j);}
    os << _T("---> Component ") << k << _T("\n");
    os << _T("mean = ") << m;
    os << _T("sigma = ")<< s;
  }
}

} // namespace STK

#endif /* STK_HDGaussianBASE_H */
