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
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_MixturePoisson_lk.h
 *  @brief In this file we implement the MixturePoisson_lk class
 **/

#ifndef STK_MIXTUREPOISSON_LK_H
#define STK_MIXTUREPOISSON_LK_H

#include "STK_MixturePoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class MixturePoisson_lk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the MixturePoisson_lk traits policy. */
template<class Array_>
struct MixtureTraits< MixturePoisson_lk<Array_> >
{
  typedef Array_                  Array;
  typedef typename Array::Type    Type;
  /** Type of the structure storing the parameters of a MixtureMixturePoisson_lk model*/
  typedef ModelParameters<Clust::Poisson_lk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixturePoisson_lk model.
 */
template<>
struct ModelParameters<Clust::Poisson_lk_>
{
    /** intensity of the variables */
    Array1D<Real> lambda_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): lambda_(nbCluster) {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : lambda_(param.lambda_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the intensity of the kth cluster and jth variable */
    inline Real const& lambda(int k, int j) const { return lambda_[k];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      for (int k = lambda_.begin(); k< lambda_.end(); ++k)
      {
        lambda_[k] = 1.;
      }
    }
};

/** @ingroup Clustering
 *  The Poisson mixture model @c MixturePoisson_lk has a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{k}} \frac{\lambda_{k}^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class MixturePoisson_lk : public MixturePoissonBase<MixturePoisson_lk<Array> >
{
  public:
    typedef MixturePoissonBase<MixturePoisson_lk<Array> > Base;
    using Base::param_;
    using Base::p_data;
    using Base::nbVariable;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    MixturePoisson_lk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    MixturePoisson_lk( MixturePoisson_lk const& model) : Base(model) {}
    /** destructor */
    ~MixturePoisson_lk() {}
    /** @return the value of lambda of the kth cluster and jth variable */
    inline Real lambda(int k, int j) const { return param_.lambda_[k];}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0., lambda = param_.lambda_[k];
      if (lambda)
      {
        for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
        { sum += Law::Poisson::lpdf(p_data()->elt(i,j), lambda);}
      }
      return sum;
    }
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** Compute the weighted probabilities. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return this->nbCluster();}
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void MixturePoisson_lk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  Real m = p_data()->template cast<Real>().mean();
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  { param_.lambda_[k] = Law::Exponential::rand(m);}
}


/* Compute the modalities probabilities */
template<class Array>
bool MixturePoisson_lk<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.lambda_[k]= (p_data()->transpose() * p_tik->col(k)).sum()
                      /(p_data()->sizeCols()*p_nk->elt(k));
  }
  return true;
}

} // namespace STK

#endif /* STK_MIXTUREPOISSON_LK_H */
