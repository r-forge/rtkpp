/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015 Serge Iovleff

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

/** @file STK_DiagGaussianBase.h
 *  @brief In this file we implement the base class for the Gaussian diagonal models
 **/

#ifndef STK_DIAGGAUSSIANBASE_H
#define STK_DIAGGAUSSIANBASE_H

#include "../STK_IMixtureModel.h"
#include "../STK_MixtureParameters.h"
#include <STatistiK/include/STK_Law_Normal.h>
#include <STatistiK/include/STK_Law_Uniform.h>

namespace STK
{


/** @ingroup Clustering
 *  Base class for the diagonal Gaussian models
 **/
template<class Derived>
class DiagGaussianBase : public IMixtureModel<Derived >
{
  public:
    typedef IMixtureModel<Derived > Base;
    typedef typename Clust::MixtureTraits<Derived>::ParamHandler ParamHandler;
    using Base::p_tik; using Base::param_;
    using Base::p_data;

  protected:
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline DiagGaussianBase( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline DiagGaussianBase( DiagGaussianBase const& model) : Base(model) {}
    /** destructor */
    inline ~DiagGaussianBase() {}

  public:
    /** @return the mean of the kth cluster and jth variable */
    inline Real mean(int k, int j) const { return this->asDerived().meanImpl(k,j);}
    /** @return the mean of the kth cluster and jth variable */
    inline Real sigma(int k, int j) const { return this->asDerived().sigmaImpl(k,j);}
    /** Initialize the parameters of the model. */
    inline void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute */
    Real impute(int i, int j) const
    {
      Real sum = 0.;
      for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
      { sum += p_tik()->elt(i,k) * mean(j,k);}
      return sum;
    }
    /** @return a simulated value for the jth variable of the ith sample
     * @param i,j indexes of the data to simulate */
    Real sample(int i, int j) const
    {
      int k = Law::Categorical::rand(p_tik()->row(i));
      return Law::Normal::rand(mean(j,k), sigma(j,k));
    }
    /** get the parameters of the model
     *  @param params the parameters of the model
     **/
    void getParameters(ArrayXX& params) const;
    /** @return the parameters of the model in an array of size (K * 2d). */
    ArrayXX getParameters() const;
    /** Write the parameters on the output stream os */
    void writeParameters(ostream& os) const
    {
      Array2DPoint<Real> m(p_data()->cols());
      Array2DPoint<Real> s(p_data()->cols());
      for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
      {
        // store sigma values in an array for a nice output
        for (int j= s.begin();  j < s.end(); ++j)
        { m[j] = mean(k,j); s[j] = sigma(k,j);}
        os << _T("---> Component ") << k << _T("\n");
        os << _T("mean = ") << m;
        os << _T("sigma = ")<< s;
      }
    }

  protected:
    PointX& mean(int k) { return param_.mean_[k];}
    /** sample randomly the mean of each component by sampling randomly a row
     *  of the data set.
     **/
    void randomMean();
    /** compute the weighted mean of a Gaussian mixture. */
    bool updateMean();
};

template<class Derived>
void DiagGaussianBase<Derived>::randomMean()
{
  // indexes array
  VectorXi indexes(p_data()->rows());
  for(int i=p_data()->beginRows(); i< p_data()->endRows(); ++i) { indexes[i] = i;}
  Range rind = p_data()->rows();
  // sample without repetition
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    // random number in [0, end-k[
    int i = (int)Law::Uniform::rand(rind.begin(), rind.end());
    // get ith individuals
    mean(k).copy(p_data()->row(indexes[i]));
    // exchange it with nth
    indexes.swap(i, rind.lastIdx());
    // decrease
    rind.decLast(1);
  }
}

template<class Derived>
bool DiagGaussianBase<Derived>::updateMean()
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
    { mean(k)[j] = p_data()->col(j).wmean(p_tik()->col(k));}
  }
  return true;
}

/* get the parameters of the model
 *  @param params the parameters of the model
 **/
template<class Derived>
void DiagGaussianBase<Derived>::getParameters(ArrayXX& params) const
{
  int nbClust = this->nbCluster();
  params.resize(2*nbClust, p_data()->cols());
  for (int k= 0; k < nbClust; ++k)
  {
    for (int j=  p_data()->beginCols();  j < p_data()->endCols(); ++j)
    {
      params(baseIdx+2*k  , j) = mean(k,j);
      params(baseIdx+2*k+1, j) = sigma(k,j);
    }
  }
}
/* @return the parameters of the model in an array of size (K * 2d). */
template<class Derived>
ArrayXX DiagGaussianBase<Derived>::getParameters() const
{
  ArrayXX params;
  int nbClust = this->nbCluster();
  params.resize(2*nbClust, p_data()->cols());
  for (int k= 0; k < nbClust; ++k)
  {
    for (int j=  p_data()->beginCols();  j < p_data()->endCols(); ++j)
    {
      params(baseIdx+2*k  , j) = mean(k,j);
      params(baseIdx+2*k+1, j) = sigma(k,j);
    }
  }
  return params;
}

} // namespace STK

#endif /* STK_DIAGGAUSSIANBASE_H */
