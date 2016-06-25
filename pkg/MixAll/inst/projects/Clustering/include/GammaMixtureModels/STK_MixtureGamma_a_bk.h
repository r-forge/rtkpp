/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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

/** @file STK_MixtureGamma_a_bk.h
 *  @brief In this file we define the Gamma_pk_a_bk and Gamma_p_a_bk models.
 **/

#ifndef STK_MIXTUREGAMMA_A_BK_H
#define STK_MIXTUREGAMMA_A_BK_H

#include "STK_MixtureGammaBase.h"
#include <STatistiK/include/STK_Law_Exponential.h>

namespace STK
{
template<class Array>class MixtureGamma_a_bk;

namespace hidden
{
/** @ingroup Clustering
 * Traits class for the MixtureGamma_a_bk traits policy
 **/
template<class Array_>
struct MixtureTraits< MixtureGamma_a_bk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a Mixture Gamma_a_bk model*/
  typedef ModelParameters<Clust::Gamma_a_bk_> Parameters;
};

} // namespace Clust

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixtureGamma_a_bk model.
 */
template<>
struct ModelParameters<Clust::Gamma_a_bk_>: public ParametersGammaBase
{
    /** shape of the variables */
    Real shape_;
    /** scales of the variables */
    Array1D<Real>  scale_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): ParametersGammaBase(nbCluster), shape_(0.), scale_(nbCluster) {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : ParametersGammaBase(param), shape_(param.shape_), scale_(param.scale_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_;}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      ParametersGammaBase::resize(range);
      shape_ = 1.;
      for (int k = scale_.begin(); k< scale_.end(); ++k)
      { scale_[k] = 1.;}
    }
};

/** @ingroup Clustering
 *  MixtureGamma_a_bk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{k}}\right)^{a-1}
 *                   \frac{e^{-x_i^j/b_{k}}}{b_{k} \, \Gamma(a)},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class MixtureGamma_a_bk : public MixtureGammaBase< MixtureGamma_a_bk<Array> >
{
  public:
    typedef MixtureGammaBase< MixtureGamma_a_bk<Array> > Base;
    using Base::param_;
    
    using Base::p_data;
    using Base::moments;
    using Base::meanjk;
    using Base::variancejk;
    using Base::meank;
    using Base::variancek;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    MixtureGamma_a_bk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    MixtureGamma_a_bk( MixtureGamma_a_bk const& model): Base(model) {}
    /** destructor */
    ~MixtureGamma_a_bk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        if (param_.shape_ && param_.scale_[k])
        { sum += Law::Gamma::lpdf(p_data()->elt(i,j), param_.shape_, param_.scale_[k]);}
      }
      return sum;
    }
    /** Initialize randomly the parameters of the Gamma mixture. The shape
     *  will be selected randomly using an exponential of parameter mean^2/variance
     *  and the scale will be selected randomly using an exponential of parameter
     *  variance/mean.
     */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** Compute the run( CArrayXX const*  p_tik, CPointX const* p_nk) . */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return this->nbCluster()+1;}
};

template<class Array>
void MixtureGamma_a_bk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  // compute moments
  this->moments(p_tik);
  Real value = 0.0;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    Real mean = meank(k), variance = variancek(k);
    // generate scales
    param_.scale_[k] = Law::Exponential::rand(variance/mean);
    value += p_nk->elt(k) * (mean*mean/variance);
  }
  param_.shape_ = STK::Law::Exponential::rand(value/this->nbSample());
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T(" Gamma_a_bk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk)  done\n");
#endif
}


/* Compute the weighted mean and the common variance. */
template<class Array>
bool MixtureGamma_a_bk<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  if (!moments(p_tik)) { return false;}
  // estimate a
  Real y =0.0, x0 = 0.0, x1 = param_.shape_;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    Real mean = meank(k);
    y  += p_nk->elt(k) * (param_.meanLog_[k] - std::log(mean)).sum();
    x0 += p_nk->elt(k) * (mean*mean/variancek(k));
  }
  y  /= (this->nbSample()*this->nbVariable());
  x0 /= this->nbSample();
  // moment estimate and oldest value
  if ((x0 <=0.) || (isNA(x0))) return false;

  // get shape
  hidden::invPsiMLog f(y);
  Real a = Algo::findZero(f, x0, x1, 1e-08);
  if (!Arithmetic<Real>::isFinite(a))
  {
#ifdef STK_MIXTURE_DEBUG
    stk_cout << "ML estimation failed in MixtureGamma_a_bjk::run( CArrayXX const*  p_tik, CPointX const* p_nk) \n";
    stk_cout << "x0 =" << x0 << _T("\n";);
    stk_cout << "f(x0) =" << f(x0) << _T("\n";);
    stk_cout << "x1 =" << x1 << _T("\n";);
    stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
    a = x0; // use moment estimate
  }
  // set values
  param_.shape_ = a;
  // estimate b
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  { param_.scale_[k] = meank(k)/a;}
  return true;
}

}  // namespace STK


#endif /* STK_MIXTUREGAMMA_A_BK_H */
