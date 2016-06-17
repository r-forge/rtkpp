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

/** @file STK_MixtureCategoricalBase.h
 *  @brief In this file we implement the base class for the Categorical
 *  diagonal models
 **/

#ifndef STK_MIXTURECATEGORICALBASE_H
#define STK_MIXTURECATEGORICALBASE_H

#include "../STK_IMixtureDensity.h"
#include <STatistiK/include/STK_Law_Categorical.h>

namespace STK
{

/** @ingroup Clustering
 *  Base class for the Categorical models
 **/
template<class Derived>
class MixtureCategoricalBase : public IMixtureDensity<Derived >
{
  protected:
    typedef IMixtureDensity<Derived> Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     *  @param nbCluster number of cluster in the model
     **/
    inline MixtureCategoricalBase( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixtureCategoricalBase( MixtureCategoricalBase const& model)
                          : Base(model), modalities_(model.modalities_) {}
    /** destructor */
    inline ~MixtureCategoricalBase() {}

  public:
    /** @return the array with the number of modalities of each columns in data set */
    inline PointXi const& nbModalities() const { return nbModalities_;}
    /** @return the range of the modalities */
    inline Range const& modalities() const { return modalities_;}
    /** @return the probability of the kth cluster, jth variable, lth modality */
    inline Real proba(int k, int j, int l) const { return param_.proba(k,j,l);}
    /** @return the probabilities of the kth cluster for the jth variable */
    inline CVectorX proba(int k, int j) const { return param_.proba(k,j);}
    /** Initialize the model. Resize the probability arrays of each component.*/
    void initializeModelImpl()
    {
      // compute the maximal number of modalities
      nbModalities_.resize(p_data()->cols());
      int amin = Arithmetic<int>::max(), amax = Arithmetic<int>::min();
      for (int j= p_data()->beginCols(); j < p_data()->endCols(); ++j)
      {
        int min = p_data()->col(j).minElt(), max = p_data()->col(j).maxElt();
        amin = std::min(amin, min); amax = std::max(amax, max);
        nbModalities_[j] = max-min+1;
      }
      // set range of the modalities
      modalities_ = _R(amin, amax);
      // resize vectors of probabilities
      param_.resize(modalities_,p_data()->cols());
    }
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param tk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    int impute(int i, int j, Weights const& tk) const;
    /** @return a simulated value for the jth variable of the ith sample
     * in the kth cluster
     * @param i,j,k indexes of the data to simulate */
    inline int rand(int i, int j, int k) const
    { return Law::Categorical::rand(proba(k,j));}

  protected:
    /** Array with the number of modalities of each columns of the data set */
    PointXi nbModalities_;
    /** range of the modalities */
    Range modalities_;
};

/* Implementation  */
template<class Derived>
template<class Weights>
int MixtureCategoricalBase<Derived>::impute(int i, int j, Weights const& tk) const
{
  int lmax = modalities_.begin();
  Real pmax = -Arithmetic<Real>::max();
  // compute for each modality the pondered probability of occurrence
  for (int l=modalities_.begin(); l< modalities_.end(); ++l)
  {
    Real p = 0.;
    for (int k= tk.begin(); k < tk.end(); ++k)
    { p += tk[k] * proba(k, j, l);}

    if (pmax < p) { pmax = p; lmax = l;}
  }
  return lmax;
}


} // namespace STK

#endif /* STK_MIXTURECATEGORICALBASE_H */
