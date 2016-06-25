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

/** @file STK_MixtureGaussian_sj.h
 *  @brief In this file we define and implement the MixtureGaussian_sj class
 **/

#ifndef STK_MIXTUREGAUSSIAN_SJ_H
#define STK_MIXTUREGAUSSIAN_SJ_H

#include "STK_MixtureDiagGaussianBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class MixtureGaussian_sj;

namespace hidden
{
/** @ingroup hidden
 *  Traits class for the Gaussian_s traits policy. */
template<class Array_>
struct MixtureTraits< MixtureGaussian_sj<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixturGaussian_sj model*/
  typedef ModelParameters<Clust::Gaussian_sj_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixturGaussian_sj model.
 */
template<>
struct ModelParameters<Clust::Gaussian_sj_>
{
    /** array of size nbCluster with the parameters mean of the variables */
    Array1D<CPointX> mean_;
    /** standard deviation of the variables */
    CPointX sigma_;
    /** default constructor
     *  @param nbCluster the number of clmass of the mixture
     **/
    ModelParameters(int nbCluster): mean_(nbCluster), sigma_() {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : mean_(param.mean_), sigma_(param.sigma_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the mean of the kth cluster and jth variable */
    inline Real const& mean(int k, int j) const { return mean_[k][j];}
    /** @return the standard deviation of the kth cluster and jth variable */
    inline Real const& sigma(int k, int j) const { return sigma_[j];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      for (int k = mean_.begin(); k< mean_.end(); ++k)
      { mean_[k].resize(range) = 0;}
      sigma_.resize(range) = 1.;
    }
};

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model MixtureGaussian_sj have a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma_j} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2\sigma_j^2}\right\}.
 * \f]
 **/
template<class Array>
class MixtureGaussian_sj : public MixtureDiagGaussianBase<MixtureGaussian_sj<Array> >
{
  public:
    typedef MixtureDiagGaussianBase<MixtureGaussian_sj<Array> > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline MixtureGaussian_sj( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixtureGaussian_sj( MixtureGaussian_sj const& model)
                      : Base(model)
    {}
    /** destructor */
    inline ~MixtureGaussian_sj() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        if (param_.sigma_[j])
        { sum += Law::Normal::lpdf(p_data()->elt(i,j), param_.mean_[k][j], param_.sigma_[j]);}
      }
      return sum;
    }
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviation
     *  will be set to 1.
     */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** Compute the weighted mean and the common standard deviation. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable()+this->nbVariable();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void MixtureGaussian_sj<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  // compute the initial mean
  this->randomMean(p_tik);
  // compute the standard deviation
  Array2DPoint<Real> variance(p_data()->cols(), 0.);
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance += p_tik->col(k).transpose()
               *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
                ).square()
                ;
  }
  // store the standard deviation
  param_.sigma_ = (variance /= this->nbSample()).sqrt();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("MixtureGaussian_sj<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk)  done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool MixtureGaussian_sj<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  // compute the means
  if (!this->updateMean(p_tik)) return false;
  // compute the standard deviation
  Array2DPoint<Real> variance(p_data()->cols(), 0.);
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance += p_tik->col(k).transpose()
               *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
                ).square()
                ;
  }
//  if (variance.nbAvailableValues() != this->nbVariable()) return false;
//  if ((variance > 0.).template cast<int>().sum() != this->nbVariable()) return false;
  // compute the standard deviation
  param_.sigma_ = (variance /= this->nbSample()).sqrt();
  return true;
}

} // namespace STK

#endif /* STK_MIXTUREGAUSSIAN_SJ_H */
