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

/** @file STK_MixtureCategorical_pjk.h
 *  @brief In this file we implement the MixtureCategorical_pjk class
 **/

#ifndef STK_MIXTURECATEGORICAL_PJK_H
#define STK_MIXTURECATEGORICAL_PJK_H

#include "STK_MixtureCategoricalBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class MixtureCategorical_pjk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the MixtureCategorical_pjk traits policy. */
template<class Array_>
struct MixtureTraits< MixtureCategorical_pjk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixtureCategorical_pjk model*/
  typedef ModelParameters<Clust::Categorical_pjk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixtureGaussian_s model.
 */
template<>
struct ModelParameters<Clust::Categorical_pjk_>
{
    /** array of size nbCluster with the probabilities of the variables in each
     *  classes */
    Array1D<CArrayXX> proba_;
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
    inline Real const& proba(int k, int j, int l) const { return proba_[k](l,j);}
    /** @return the probabilities of the kth cluster for the jth variable */
    inline CVectorX proba(int k, int j) const { return proba_[k].col(j);}

    /** resize the set of parameter */
    inline void resize(Range const& rangeModalities, Range const& rangeCols)
    {
      for (int k = proba_.begin(); k< proba_.end(); ++k)
      { proba_[k].resize(rangeModalities, rangeCols) = 1./rangeModalities.size();}
    }
};

/** @ingroup Clustering
 *  The diagonal Categorical mixture model @c MixtureCategorical_pjk is
 *  the most general diagonal Categorical model and have a probability
 *  function of the form
 * \f[
 *    P(\mathbf{x}=(l_1,\ldots,l_d)|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d p_{kl_j}^j.
 * \f]
 **/
template<class Array>
class MixtureCategorical_pjk : public MixtureCategoricalBase<MixtureCategorical_pjk<Array> >
{
  public:
    typedef MixtureCategoricalBase<MixtureCategorical_pjk<Array> > Base;
    using Base::param_;
    using Base::p_data;
    using Base::modalities_;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    MixtureCategorical_pjk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    MixtureCategorical_pjk( MixtureCategorical_pjk const& model): Base(model) {}
    /** destructor */
    ~MixtureCategorical_pjk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      { // what to do if the probability is zero but a sample get this modality
        // for now, just ignore it (it's possible if tik_(i,k) == 0)
        Real prob= param_.proba_[k](p_data()->elt(i,j), j);
        if (prob) { sum += std::log(prob);}
       }
      return sum;
    }
    /** Initialize randomly the parameters of the Categorical mixture. */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** Compute the weighted probabilities. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*((this->nbModalities_-1).sum());}
};

/* Initialize randomly the parameters of the Categorical mixture. */
template<class Array>
void MixtureCategorical_pjk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  for (int k = p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.proba_[k].randUnif();
    for (int j=param_.proba_[k].beginCols(); j< param_.proba_[k].endCols(); ++j)
    { param_.proba_[k].col(j) /= param_.proba_[k].col(j).sum();}
  }
}


/* Compute the modalities probabilities */
template<class Array>
bool MixtureCategorical_pjk<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk)
{
  for (int k = p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.proba_[k] = 0.;
    for (int j = p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      // count the number of modalities weighted by the tik
      for (int i = p_data()->beginRows(); i < p_data()->endRows(); ++i)
      { param_.proba_[k](p_data()->elt(i, j), j) += p_tik->elt(i, k);}
      // normalize the probabilities
      Real sum = param_.proba_[k].col(j).sum();
      if (sum) { param_.proba_[k].col(j) /= sum;}
    }
  }
  return true;
}

} // namespace STK

#endif /* STK_MIXTURECATEGORICAL_PJK_H */
