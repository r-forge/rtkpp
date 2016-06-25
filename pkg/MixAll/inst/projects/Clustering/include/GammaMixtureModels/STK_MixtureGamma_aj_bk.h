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

/** @file STK_MixtureGamma_aj_bk.h
 *  @brief In this file we define the Gamma_pk_aj_bk and Gamma_p_aj_bk models.
 **/

#ifndef STK_MIXTUREGAMMA_AJ_BK_H
#define STK_MIXTUREGAMMA_AJ_BK_H

#include "STK_MixtureGammaBase.h"
#include <STatistiK/include/STK_Law_Exponential.h>

#define MAXITER 400
#define TOL 1e-8

namespace STK
{
template<class Array>class MixtureGamma_aj_bk;

namespace hidden
{
/** @ingroup Clustering
 * Traits class for the MixtureGamma_aj_bk traits policy
 **/
template<class Array_>
struct MixtureTraits< MixtureGamma_aj_bk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixtureGamma_aj_bk model*/
  typedef ModelParameters<Clust::Gamma_aj_bk_> Parameters;
};

} // namespace Clust

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixtureGamma_aj_bk_ model.
 */
template<>
struct ModelParameters<Clust::Gamma_aj_bk_>: public ParametersGammaBase
{
    /** shapes of the variables */
    CPointX shape_;
    /** scales of the variables */
    Array1D<Real>  scale_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): ParametersGammaBase(nbCluster), shape_(), scale_(nbCluster) {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : ParametersGammaBase(param), shape_(param.shape_), scale_(param.scale_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the mean of the kth cluster and jth variable */
    inline Real const& shape(int k, int j) const { return shape_[j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& scale(int k, int j) const { return scale_[k];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      ParametersGammaBase::resize(range);
      shape_.resize(range) = 1.;
      for (int k = scale_.begin(); k< scale_.end(); ++k)
      { scale_[k] = 1.;}
    }
};

/** @ingroup Clustering
 *  MixtureGamma_aj_bk is a mixture model of the following form
 * \f[
 *     f(\mathbf{x}_i|\theta) = \sum_{k=1}^K p_k
 *     \prod_{j=1}^p\left(\frac{x_i^j}{b_{k}}\right)^{a_{j}-1}
 *                   \frac{e^{-x_i^j/b_{k}}}{b_{k} \, \Gamma(a_{j})},
 *      \quad x_i^j>0, \quad i=1,\ldots,n.
 * \f]
 **/
template<class Array>
class MixtureGamma_aj_bk : public MixtureGammaBase< MixtureGamma_aj_bk<Array> >
{
  public:
    typedef MixtureGammaBase< MixtureGamma_aj_bk<Array> > Base;
    using Base::param_;
    
    using Base::p_data;
    using Base::meanjk;
    using Base::variancejk;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    MixtureGamma_aj_bk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    MixtureGamma_aj_bk( MixtureGamma_aj_bk const& model): Base(model) {}
    /** destructor */
    ~MixtureGamma_aj_bk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        if (param_.shape_[j] && param_.scale_[k])
        { sum += Law::Gamma::lpdf(p_data()->elt(i,j), param_.shape_[j], param_.scale_[k]);}
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
    inline int computeNbFreeParameters() const
    { return(this->nbCluster()+this->nbVariable());}
};

/* Initialize randomly the parameters of the gamma mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void MixtureGamma_aj_bk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
    // compute moments
    this->moments(p_tik);
  // simulate aj
  for (int j=p_data()->beginCols(); j < p_data()->endCols(); ++j)
  {
    Real value= 0.;
    for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
    {
      Real mean = meanjk(j,k), variance = variancejk(j,k);
      value += p_nk->elt(k) * mean*mean/variance;
    }
    param_.shape_[j] = Law::Exponential::rand(value/(this->nbSample()));
  }
  // simulates bk
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  { param_.scale_[k] = Law::Exponential::rand((this->variancek(k)/this->meank(k)));}
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T(" Gamma_aj_bk<Array>::randomInit done\n");
#endif
}

/* Compute the weighted mean and the common variance. */
template<class Array>
bool MixtureGamma_aj_bk<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  if (!this->moments(p_tik)) { return false;}
  // start estimations of the ajk and bj
  Real qvalue = this->qValue(p_tik, p_nk);
  // enter iterative algorithm
  int iter;
  for(iter = 0; iter<MAXITER; ++iter)
  {
    // compute aj
    for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
    {
      // moment estimate and oldest value
      Real y = 0, x0 = 0;
      for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
      {
        Real mean = meanjk(j,k), variance = variancejk(j,k);
        y  += p_nk->elt(k) * (param_.meanLog_[k][j] - std::log(param_.scale_[k]));
        x0 += p_nk->elt(k) * mean*mean/variance;
      }
      y /= this->nbSample();
      x0/= this->nbSample();
      Real x1 = param_.shape_[j];
      if ((x0 <=0.) || !Arithmetic<Real>::isFinite(x0)) return false;
      // compute shape
      hidden::invPsi f(y);
      Real a =  Algo::findZero(f, x0, x1, TOL);

      if (!Arithmetic<Real>::isFinite(a))
      {
        param_.shape_[j] = x0; // use moment estimate
#ifdef STK_MIXTURE_DEBUG
        stk_cout << _T("ML estimation failed in MixtureGamma_ajk_bj::run( CArrayXX const*  p_tik, CPointX const* p_nk) \n");
        stk_cout << "x0 =" << x0 << _T("\n";);
        stk_cout << "f(x0) =" << f(x0) << _T("\n";);
        stk_cout << "x1 =" << x1 << _T("\n";);
        stk_cout << "f(x1) =" << f(x1) << _T("\n";);
#endif
      }
      else { param_.shape_[j] = a;}
      // compute bk
      Real sum = param_.shape_.sum();
      for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
      { // update bk
        param_.scale_[k] = param_.mean_[k].sum()/sum;
      }
    }
  // check convergence
    Real value = this->qValue(p_tik, p_nk);
#ifdef STK_MIXTURE_VERBOSE
  if (value < qvalue)
  {
    stk_cout << _T("In MixtureGamma_aj_bk::run( CArrayXX const*  p_tik, CPointX const* p_nk) : run( CArrayXX const*  p_tik, CPointX const* p_nk)  diverge\n");
    stk_cout << _T("New value =") << value << _T(", qvalue =") << qvalue << _T("\n");
  }
#endif
    if ((value - qvalue) < TOL) break;
    qvalue = value;
  }
#ifdef STK_MIXTURE_VERBOSE
  if (iter == MAXITER)
  {
    stk_cout << _T("In MixtureGamma_aj_bk::run( CArrayXX const*  p_tik, CPointX const* p_nk) : run( CArrayXX const*  p_tik, CPointX const* p_nk)  did not converge\n");
    stk_cout << _T("qvalue =") << qvalue << _T("\n");
  }
#endif
  return true;
}

}  // namespace STK

#undef MAXITER
#undef TOL

#endif /* STK_MIXTUREGAMMA_AJ_BK_H */
