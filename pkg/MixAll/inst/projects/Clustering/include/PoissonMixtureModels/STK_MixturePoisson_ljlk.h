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

/** @file STK_MixturePoisson_ljlk.h
 *  @brief In this file we implement the MixturePoisson_ljlk class
 **/

#ifndef STK_MIXTUREPOISSON_LJLK_H
#define STK_MIXTUREPOISSON_LJLK_H

#include "STK_MixturePoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class MixturePoisson_ljlk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the MixturePoisson_ljlk traits policy. */
template<class Array_>
struct MixtureTraits< MixturePoisson_ljlk<Array_> >
{
  typedef Array_ Array;
  typedef typename Array::Type    Type;
  /** Type of the structure storing the parameters of a MixtureMixturePoisson_ljlk model*/
  typedef ModelParameters<Clust::Poisson_ljlk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixturePoisson_ljlk model.
 */
template<>
struct ModelParameters<Clust::Poisson_ljlk_>
{
    /** intensity of the variables by class */
    Array1D<Real> lambdak_;
    /** intensity of the variables by variables */
    CVectorX lambdaj_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): lambdak_(nbCluster), lambdaj_() {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : lambdak_(param.lambdak_), lambdaj_(param.lambdaj_)
    {}
    /** @return the intensity of the kth cluster and jth variable */
    inline Real lambda(int k, int j) const { return lambdak_[k] * lambdaj_[j];}
    /** destructor */
    ~ModelParameters() {}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      for (int k = lambdak_.begin(); k< lambdak_.end(); ++k)
      {
        lambdak_[k] = 1.;
      }
      lambdaj_.resize(range) = 1;
    }
};

/** @ingroup Clustering
 *  The Poisson mixture model @c MixturePoisson_ljlk is a Poisson model
 *  with a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{j}\lambda_{k}} \frac{(\lambda_{j}\lambda_{k})^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class MixturePoisson_ljlk : public MixturePoissonBase<MixturePoisson_ljlk<Array> >
{
  public:
    typedef MixturePoissonBase<MixturePoisson_ljlk<Array> > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    MixturePoisson_ljlk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    MixturePoisson_ljlk( MixturePoisson_ljlk const& model): Base(model) {}
    /** destructor */
    ~MixturePoisson_ljlk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        Real lambda = param_.lambdak_[k]*param_.lambdaj_[j];
        if (lambda)
        { sum += Law::Poisson::lpdf(p_data()->elt(i,j), lambda);}
      }
      return sum;
    }
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** Compute the weighted probabilities. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()+this->nbVariable();}
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void MixturePoisson_ljlk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
  {
    Real m = p_data()->col(j).template cast<Real>().mean();
    for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
    {
      param_.lambdak_[k] = Law::Exponential::rand(m)/param_.lambdaj_[j];
    }
  }
}


/* Compute the lambdas */
template<class Array>
bool MixturePoisson_ljlk<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  param_.lambdaj_ = (Stat::sumByRow(*p_tik).transpose() * (*p_data()))
                     /(Stat::sumByRow(*p_tik) * Stat::sumByRow(*p_data())).sum();
  param_.lambdak_ = Stat::sumByRow(*p_data()).transpose() * (*p_tik)/(*p_nk);
  return true;
}

} // namespace STK

#endif /* STK_MIXTUREPOISSON_LJLK_H */
