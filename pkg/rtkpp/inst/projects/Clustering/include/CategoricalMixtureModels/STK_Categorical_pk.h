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

/** @file STK_Categorical_pk.h
 *  @brief In this file we define the Categorical_pk model
 **/

#ifndef STK_CATEGORICAL_PK_H
#define STK_CATEGORICAL_PK_H

#include "STK_CategoricalBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class Categorical_pk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the Categorical_pk traits policy. */
template<class _Array>
struct MixtureTraits< Categorical_pk<_Array> >
{
  typedef _Array Array;
  typedef typename Array::Type Type;
  typedef Array2D<Real>   Param;
  typedef ParametersHandler<Clust::Categorical_pk_> ParamHandler;
};

} // namespace hidden

/** Specialization of the ParametersHandler struct for Categorical_pk model */
template <>
struct ParametersHandler<Clust::Categorical_pk_>
{
  /** default constructor */
  ParametersHandler(int nbCluster):proba_(nbCluster) {}
  /** copy constructor */
  ParametersHandler(ParametersHandler const& model):proba_(model.proba_) {}
  /** Initialize the parameters of the model.
   *  This function initialize the parameter lambda and the statistics.
   **/
  void resize( Range const& rangeModalities, Range const& rangeData)
  {
    proba_.resize(rangeModalities);
    proba_.initialize(1./rangeModalities.size());
  }
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { proba_.storeIntermediateResults(iteration);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { proba_.releaseIntermediateResults();}
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters() { proba_.setParameters();}
  /** Vector and statistics of the probabilities */
  MixtureParametersSet<VectorX> proba_;
};

/** @ingroup Clustering
 *  The diagonal Categorical mixture model @c Categorical_pk is
 *  the most general diagonal Categorical model and have a density function of the
 *  form
 * \f[
 *  P(\mathbf{x}=(l_1,\ldots,l_d)|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d p_{kl_j}.
 * \f]
 **/
template<class Array>
class Categorical_pk : public CategoricalBase<Categorical_pk<Array> >
{
  public:
    typedef CategoricalBase<Categorical_pk<Array> > Base;
    using Base::p_tik; using Base::param_;
    using Base::p_data;
    using Base::modalities_;
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Categorical_pk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Categorical_pk( Categorical_pk const& model): Base(model) {}
    /** destructor */
    inline ~Categorical_pk() {}
    /** @return the probability of the kth cluster, jth variable, lth modality */
    inline Real probaImpl(int k, int j, int l) const { return param_.proba_[k][l];}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0., prob;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        if ( (prob = param_.proba_[k][p_data()->elt(i,j)]) <= 0.) return -Arithmetic<Real>::infinity();
        sum += std::log(prob);
       }
      return sum;
    }
    /** Initialize randomly the parameters of the Categorical mixture.
     *  Probabilities will be choosen uniformly.
     */
    void randomInit();
    /** Compute the weighted proportions of each class. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*(this->modalities_.size()-1);}
};

/* Initialize randomly the parameters of the Categorical mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Categorical_pk<Array>::randomInit()
{
  for (int k = p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.proba_[k].randUnif();
    param_.proba_[k] /= param_.proba_[k].sum();
  }
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Categorical_pk<Array>::mStep()
{
  for (int k = p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.proba_[k] = 0.;
    for (int j = p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      for (int i = p_tik()->beginRows(); i < p_tik()->endRows(); ++i)
      { param_.proba_[k][(*p_data())(i, j)] += (*p_tik())(i, k);}
    }
    Real sum = param_.proba_[k].sum();
    if (sum<=0.) return false;
    param_.proba_[k] /= sum;
  }
  return true;
}

} // namespace STK

#endif /* STK_CATEGORICAL_PK_H */
