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

/* Project:  stkpp::Clustering
 * created on: Dec 9, 2014
 * Authors: Serge Iovleff
 **/

/** @file STK_MixturePoissonBase.h
 *  @brief In this file we implement the base class for the poisson mixture models
 **/

#ifndef STK_MIXTUREPOISSONBASE_H
#define STK_MIXTUREPOISSONBASE_H

#include "../STK_IMixtureDensity.h"
#include <STatistiK/include/STK_Law_Poisson.h>

namespace STK
{
/** @ingroup Clustering
 *  Base class for the Poisson models
 **/
template<class Derived>
class MixturePoissonBase: public IMixtureDensity<Derived >
{
  public:
    typedef IMixtureDensity<Derived > Base;
    using Base::param_;
    using Base::p_data;
    using Base::nbCluster;

  protected:
    /** default constructor
     *  @param nbCluster number of cluster in the model
     **/
    inline MixturePoissonBase( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixturePoissonBase( MixturePoissonBase const& model) : Base(model) {}
    /** destructor */
    inline ~MixturePoissonBase() {}

  public:
    /** @return the value of lambda of the kth cluster and jth variable */
    inline Real lambda(int k, int j) const { return param_.lambda(k,j);}
    /** Initialize the parameters of the model. */
    void initializeModelImpl() { param_.resize(p_data()->cols());}
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param pk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    int impute(int i, int j, Weights const& pk) const;
    /** @return a value to impute for the jth variable of the ith sample*/
    Real impute(int i, int j, CArrayXX const*  p_tik) const
    {
      Real sum = 0.;
      for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
      { sum += p_tik->elt(i,k) * lambda(k,j);}
      return sum;
    }
    /** @return a simulated value for the jth variable of the ith sample
     *  in the kth cluster.
     *  @param i,j,k indexes of the data to simulate
     **/
    inline int rand(int i, int j, int k) const
    { return Law::Poisson::rand(lambda(k,j));}
};

/* Implementation  */
template<class Derived>
template<class Weights>
int MixturePoissonBase<Derived>::impute(int i, int j, Weights const& pk) const
{
  Real sum = 0.;
  for (int k= pk.begin(); k < pk.end(); ++k)
  { sum += pk[k] * lambda(k,j);}
  return std::floor(sum+0.5);
}


} // namespace STK

#endif /* STK_MIXTUREPOISSONBASE_H */
