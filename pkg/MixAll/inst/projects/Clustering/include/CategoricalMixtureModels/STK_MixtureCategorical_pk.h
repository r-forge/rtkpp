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

/** @file STK_MixtureCategorical_pk.h
 *  @brief In this file we define the MixtureCategorical_pk model
 **/

#ifndef STK_CATEGORICAL_PK_H
#define STK_CATEGORICAL_PK_H

#include "STK_MixtureCategoricalBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class MixtureCategorical_pk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the MixtureCategorical_pk traits policy. */
template<class Array_>
struct MixtureTraits< MixtureCategorical_pk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixtureCategorical_pjk model*/
  typedef ModelParameters<Clust::Categorical_pk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixtureCategorical_pk model.
 */
template<>
struct ModelParameters<Clust::Categorical_pk_>
{
    /** array of size nbCluster with the probabilities of the variables */
    Array1D<CVectorX> proba_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): proba_(nbCluster) {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : proba_(param.proba_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the probability of the kth cluster, jth variable, lth modality */
    inline Real const& proba(int k, int j, int l) const { return proba_[k][l];}
    /** @return the probabilities of the kth cluster for the jth variable */
    inline CVectorX proba(int k, int j) const { return proba_[k];}
    /** resize the set of parameter */
    inline void resize(Range const& rangeModalities, Range const& range)
    {
      for (int k = proba_.begin(); k< proba_.end(); ++k)
      { proba_[k].resize(rangeModalities) = 1./rangeModalities.size();}
    }
};

/** @ingroup Clustering
 *  The diagonal Categorical mixture model @c MixtureCategorical_pk is
 *  the most general diagonal Categorical model and have a density function of the
 *  form
 * \f[
 *  P(\mathbf{x}=(l_1,\ldots,l_d)|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d p_{kl_j}.
 * \f]
 **/
template<class Array>
class MixtureCategorical_pk : public MixtureCategoricalBase<MixtureCategorical_pk<Array> >
{
  public:
    typedef MixtureCategoricalBase<MixtureCategorical_pk<Array> > Base;
    using Base::param_;
    using Base::p_data;
    using Base::modalities_;
    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline MixtureCategorical_pk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixtureCategorical_pk( MixtureCategorical_pk const& model): Base(model) {}
    /** destructor */
    inline ~MixtureCategorical_pk() {}
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
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** Compute the weighted proportions of each class. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*(this->modalities_.size()-1);}
};

/* Initialize randomly the parameters of the Categorical mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void MixtureCategorical_pk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  for (int k = p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.proba_[k].randUnif();
    param_.proba_[k] /= param_.proba_[k].sum();
  }
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool MixtureCategorical_pk<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  for (int k = p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.proba_[k] = 0.;
    for (int j = p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      for (int i = p_tik->beginRows(); i < p_tik->endRows(); ++i)
      { param_.proba_[k][(*p_data())(i, j)] += (*p_tik)(i, k);}
    }
    Real sum = param_.proba_[k].sum();
    if (sum<=0.) return false;
    param_.proba_[k] /= sum;
  }
  return true;
}

} // namespace STK

#endif /* STK_CATEGORICAL_PK_H */
