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
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_Poisson_ljlk.h
 *  @brief In this file we implement the Poisson_ljlk class
 **/

#ifndef STK_POISSON_LJLK_H
#define STK_POISSON_LJLK_H

#include "STK_PoissonBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Poisson_ljlk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Poisson_ljlk traits policy. */
template<class _Array>
struct MixtureTraits< Poisson_ljlk<_Array> >
{
  typedef _Array                  Array;
  typedef typename Array::Type    Type;
  typedef Array2D<Real>           Param;
  typedef ParametersHandler<Clust::Poisson_ljlk_> ParamHandler;
};

} // namespace hidden

/** Specialization of the ParametersHandler struct for Poisson_ljlk model */
template <>
struct ParametersHandler<Clust::Poisson_ljlk_>
{
  /** Array of the class rates */
  MixtureParametersSet<Real> lambdak_;
  /** Array of the variables rates */
  MixtureParameters<PointX> lambdaj_;
  /** default constructor. All lambdas are initialized to 1. */
  inline ParametersHandler( int nbCluster)
                          : lambdak_(nbCluster), lambdaj_()
  {}
  /** copy constructor.
   * @param param the parameters to copy.
   **/
  inline ParametersHandler( ParametersHandler const& param)
                          : lambdak_(param.lambdak_), lambdaj_(param.lambdaj_)
  {}
  /** destructor */
  inline ~ParametersHandler() {}
  /** Initialize the parameters of the model. */
  inline void resize(Range const& range)
  {
    lambdak_.initialize(1.);
    lambdaj_.resize(range);
    lambdaj_.initialize(1.);
  }
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { lambdak_.storeIntermediateResults(iteration); lambdaj_.storeIntermediateResults(iteration);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { lambdak_.releaseIntermediateResults(); lambdaj_.releaseIntermediateResults();}
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters() { lambdak_.setParameters(); lambdaj_.setParameters();}
};

/** @ingroup Clustering
 *  The Poisson mixture model @c Poisson_ljlk is a Poisson model
 *  with a probability function of the form
 * \f[
 *    P(\mathbf{x}=(n_1,\ldots,n_d)|\theta)
 *     = \sum_{k=1}^K p_k \prod_{j=1}^d e^{-\lambda_{j}\lambda_{k}} \frac{(\lambda_{j}\lambda_{k})^{n_j}}{n_j!}.
 * \f]
 **/
template<class Array>
class Poisson_ljlk : public PoissonBase<Poisson_ljlk<Array> >
{
  public:
    typedef PoissonBase<Poisson_ljlk<Array> > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_nk;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    Poisson_ljlk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    Poisson_ljlk( Poisson_ljlk const& model)
                : Base(model) {}
    /** destructor */
    inline ~Poisson_ljlk() {}
    /** @return the value of lambda of the kth cluster and jth variable */
    /** @return the value of lambda of the kth cluster and jth variable */
    inline Real lambdaImpl(int k, int j) const { return param_.lambdak_[k] * param_.lambdaj_()[j];}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      { sum += Law::Poisson::lpdf(p_data()->elt(i,j), param_.lambdak_[k]*param_.lambdaj_()[j]);}
      return sum;
    }
    /** Initialize randomly the parameters of the Poisson mixture. */
    void randomInit();
    /** Compute the weighted probabilities. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()+this->nbVariable();}
};

/* Initialize randomly the parameters of the Poisson mixture. */
template<class Array>
void Poisson_ljlk<Array>::randomInit()
{
  for (int j=p_data()->beginCols(); j< p_data()->endCols(); ++j)
  {
    Real m = p_data()->col(j).template cast<Real>().mean();
    for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
    {
      param_.lambdak_[k] = Law::Exponential::rand(m)/param_.lambdaj_()[j];
    }
  }
}


/* Compute the lambdas */
template<class Array>
bool Poisson_ljlk<Array>::mStep()
{
  param_.lambdaj_() = (Stat::sumByRow(*p_tik()).transpose() * (*p_data()))
                     /(Stat::sumByRow(*p_tik()) * Stat::sumByRow(*p_data())).sum();
  param_.lambdak_() = Stat::sumByRow(*p_data()).transpose() * (*p_tik())/(*p_nk());
  return true;
}

} // namespace STK

#endif /* STK_POISSON_LJLK_H */
