/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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
 * Project: stkpp::Clustering
 * created on: 5 sept. 2013
 * Author:  iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_Gamma_aj_bk.h
 *  @brief In this file we define the Gamma_pk_aj_bk and Gamma_p_aj_bk models.
 **/

#ifndef STK_GAMMA_AJ_BK_H
#define STK_GAMMA_AJ_BK_H

#include "STK_GammaBase.h"

#include "../../../STatistiK/include/STK_Law_Exponential.h"

namespace STK
{
template<class Array>class Gamma_aj_bk;

namespace Clust
{
/** @ingroup Clustering
 * Traits class for the Gamma_aj_bk traits policy
 **/
template<class _Array>
struct MixtureTraits< Gamma_aj_bk<_Array> >
{
  typedef _Array Array;
  typedef typename Array::Type Type;
  typedef MixtureComponent<_Array, Gamma_aj_bk_Parameters> Component;
  typedef Gamma_aj_bk_Parameters        Parameters;
  typedef Array2D<Real>        Param;
};

} // namespace Clust

/** @ingroup Clustering
 *  Gamma_aj_bk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{k}}\right)^{a_{j}-1}
 *                   \frac{e^{-x_i^j/b_{k}}}{b_{k} \, \Gamma(a_{j})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class Gamma_aj_bk : public GammaBase< Gamma_aj_bk<Array> >
{
  public:
    typedef typename Clust::MixtureTraits< Gamma_aj_bk<Array> >::Component Component;
    typedef typename Clust::MixtureTraits< Gamma_aj_bk<Array> >::Parameters Parameters;
    typedef GammaBase< Gamma_aj_bk<Array> > Base;

    using Base::p_tik;
    using Base::p_data;
    using Base::p_param;
    using Base::components;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline Gamma_aj_bk( int nbCluster) : Base(nbCluster), shape_() {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline Gamma_aj_bk( Gamma_aj_bk const& model) : Base(model), shape_(model.shape_) {}
    /** destructor */
    inline ~Gamma_aj_bk() {}
    /** Initialize the component of the model.
     *  This function have to be called prior to any used of the class.
     *  In this interface, the @c initializeModel() method call the base
     *  class IMixtureModel::initializeModel() and for all the
     *  components create the parameters.
     **/
    void initializeModel()
    {
      Base::initializeModel();
      shape_.resize(p_data()->cols());
      shape_ = 1.;
      for (int k= baseIdx; k < components().end(); ++k)
      { p_param(k)->p_shape_ = &shape_;}
    }
    /** initialize shape and scale parameters using weighted moment estimates.*/
    inline bool initializeStep() { return mStep();};
    /** Initialize randomly the parameters of the Gamma mixture. The shape
     *  will be selected randomly using an exponential of parameter mean^2/variance
     *  and the scale will be selected randomly using an exponential of parameter
     *  variance/mean.
     */
    void randomInit();
    /** Compute the mStep. */
    bool mStep();
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return(this->nbCluster()+this->nbVariable());}

  protected:
    /** common scale */
    Array2DPoint<Real> shape_;
};

/* Initialize randomly the parameters of the gamma mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void Gamma_aj_bk<Array>::randomInit()
{
    // simulates ak
    for (int k= baseIdx; k < components().end(); ++k)
    {
      Real value= 0.;
      for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
      {
        Real mean = meanjk(j,k), variance = variancejk(j,k);
        value += variance/mean;
      }
      p_param(k)->scale_ = Law::Exponential::rand(value/this->nbVariable());
    }
    // simulate bj
    for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
    {
      Real value= 0.;
      for (int k= baseIdx; k < components().end(); ++k)
      {
        Real mean = meanjk(j,k), variance = variancejk(j,k);
        value += p_param(k)->tk_ * mean*mean/variance;
      }
      shape_[j] = Law::Exponential::rand(value/this->nbSample());
    }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Gamma_aj_bk<Array>::randomInit done\n");
  this->writeParameters(stk_cout);
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool Gamma_aj_bk<Array>::mStep()
{
  if (!this->moments()) { return false;}
  return true;
}

}  // namespace STK


#endif /* STK_GAMMA_AJ_BK_H */
