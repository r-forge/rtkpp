/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IMixtureLearner.h
 *  @brief In this file we define the interface base class for learners.
 **/

#ifndef STK_IMIXTURELEARNER_H
#define STK_IMIXTURELEARNER_H

#include "STK_IMixtureStatModel.h"
#include <STatistiK/include/STK_Stat_Functors.h> // for sumByCol

namespace STK
{

/** @ingroup Clustering
 *  @brief Base class for Learner of a Mixture mixed model.
 *
 * In this interface we assume there is an underline generative model that will
 * be estimated using a MCMC algorithm.
 *
 * The pure virtual function to implement in derived class are
 * @code
 *   virtual void paramUpdateStep() = 0;
 * @endcode
 *
 * @sa IMixtureComposer
 */
class IMixtureLearner: public IMixtureStatModel //, public IMixture
{
  protected:
    /** Constructor.
     *  @param nbCluster,nbSample number of clusters and samples
     **/
    IMixtureLearner( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to clone
     **/
    IMixtureLearner( IMixtureLearner const& model);

  public:
    /** destructor */
    virtual ~IMixtureLearner();

    /** @return the state of the model*/
    inline Clust::modelState state() const { return state_;}
    /** set the state of the model : should be used by any strategy*/
    inline void setState(Clust::modelState state) { state_ = state;}

    // pure virtual
    /** Compute the model parameters given the current
     *  mixture parameters and imputation/simulation of the missing values.
     **/
    virtual void paramUpdateStep() = 0;

    // not virtual
    /** set the mixture parameters using the posterior probabilities.
     *  Proportions, numbers in each class and class labels are computed
     *  using these probabilities.
     **/
    template<class Array>
    void setMixtureParameters( Array const& tik);
    /** set the mixture parameters giving the posterior probabilities and
     *  the proportions.
     *  Numbers in each class and class labels are computed using the
     *  posterior probabilities.
     **/
    template<class Array, class RowVector>
    void setMixtureParameters( Array const& tik, RowVector const& pk);
    /** set the mixture parameters using the given class labels.
     *  Posterior probabilities, numbers in each class and proportions are
     *  computed using the class labels.
     **/
    template<class ColVector>
    void setClassLabels( ColVector const& zi);
    /** set the mixture parameters using the class labels and giving the
     *  proportions.
     **/
    template<class ColVector, class RowVector>
    void setClassLabels( ColVector const& zi, RowVector const& pk);

  private:
    /** state of the model*/
    Clust::modelState state_;
};

/* set the mixture parameters using the posterior probabilities.
 **/
template<class Array>
void IMixtureLearner::setMixtureParameters(Array const& tik)
{
  setNbSample(tik_.sizeRows());
  setNbCluster(tik_.sizeCols());
  tik_ = tik;
  nk_ = Stat::sumByCol(tik_);
  pk_ = nk_ / nbSample();
  for (int i = tik_.beginCols(); i< tik_.endCols(); ++i)
  {
    int k;
    tik_.row(i).maxElt(k);
    zi_[i] = k;
  }
}

/* set the mixture parameters using the posterior probabilities.
 **/
template<class Array, class RowVector>
void IMixtureLearner::setMixtureParameters(Array const& tik, RowVector const& pk)
{
#ifdef STK_MIXTURE_DEBUG
  if (pk_.size() != tik_.sizeCols())
  { STKRUNTIME_ERROR_2ARG(IMixtureLearner::setMixtureParameters,pk_.size(),tik_.sizeCols(),Numbers of class in tik and pk differ);}
#endif
  setNbSample(tik_.sizeRows());
  setNbCluster(tik_.sizeCols());
  tik_ = tik;
  pk_  = pk;
  nk_  = Stat::sumByCol(tik_);
  for (int i = tik_.beginCols(); i< tik_.endCols(); ++i)
  {
    int k;
    tik_.row(i).maxElt(k);
    zi_[i] = k;
  }
}

/* set the mixture parameters using the given class labels.
 *  Posterior probabilities, numbers in each class  are computed using these
 *  class labels.
 **/
template<class ColVector>
void IMixtureLearner::setClassLabels( ColVector const& zi)
{
  zi_  = zi;
  tik_ = 0.;
  for (int i=tik_.beginRows(); i < tik_.endRows(); i++)
  { tik_(i, zi_[i]) = 1.;}
  // count the number of individuals in each class
  nk_ = Stat::sumByCol(tik_);
  pk_ = nk_/nbSample();
}


/* set the mixture parameters using the class labels and giving the
 *  proportions.
 **/
template<class ColVector, class RowVector>
void IMixtureLearner::setClassLabels( ColVector const& zi, RowVector const& pk)
{
  pk_  = pk;
  zi_  = zi;
  tik_ = 0.;
  for (int i=tik_.beginRows(); i < tik_.endRows(); i++)
  { tik_(i, zi_[i]) = 1.;}
  // count the number of individuals in each class
  nk_ = Stat::sumByCol(tik_);
}

} // namespace STK

#endif /* STK_IMIXTURELEARNER_H */

