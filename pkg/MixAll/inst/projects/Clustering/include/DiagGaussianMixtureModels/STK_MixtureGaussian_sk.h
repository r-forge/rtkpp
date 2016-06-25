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

/** @file STK_MixtureGaussian_sk.h
 *  @brief In this file we define the MixtureGaussian_sk model
 **/

#ifndef STK_MIXTUREGAUSSIAN_SK_H
#define STK_MIXTUREGAUSSIAN_SK_H

#include "STK_MixtureDiagGaussianBase.h"

namespace STK
{

//forward declaration, to allow for recursive template
template<class Array>class MixtureGaussian_sk;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the MixtureGaussian_sk traits policy. */
template<class Array_>
struct MixtureTraits< MixtureGaussian_sk<Array_> >
{
  typedef Array_ Array;
  /** Type of the structure storing the parameters of a MixturGaussian_sk model*/
  typedef ModelParameters<Clust::Gaussian_sk_> Parameters;
};

} // namespace hidden

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixturGaussian_sj model.
 */
template<>
struct ModelParameters<Clust::Gaussian_sk_>
{
    /** array of size nbCluster with the parameters mean of the variables */
    Array1D<CPointX> mean_;
    /** standard deviation of the variables */
    Array1D<Real> sigma_;
    /** default constructor
     *  @param nbCluster the number of clmass of the mixture
     **/
    ModelParameters(int nbCluster): mean_(nbCluster), sigma_(nbCluster) {}
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
    inline Real const& sigma(int k, int j) const { return sigma_[k];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      for (int k = mean_.begin(); k< mean_.end(); ++k)
      { mean_[k].resize(range) =0; sigma_[k] = 1.;}
    }
};

/** @ingroup Clustering
 *  The diagonal Gaussian mixture model @c MixtureGaussian_sk assumes an equal standard
 *  deviation in each cluster and has a density function of the form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma_{k}} \exp\left\{-\frac{(x^j-\mu^j_k)^2}{2(\sigma_{k})^2}\right\}.
 * \f]
 **/
template<class Array>
class MixtureGaussian_sk : public MixtureDiagGaussianBase<MixtureGaussian_sk<Array> >
{
  public:
    typedef MixtureDiagGaussianBase<MixtureGaussian_sk<Array> > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline MixtureGaussian_sk( int nbCluster) : Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixtureGaussian_sk( MixtureGaussian_sk const& model) : Base(model) {}
    /** destructor */
    inline ~MixtureGaussian_sk() {}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      Real sum =0.;
      for (int j=p_data()->beginCols(); j<p_data()->endCols(); ++j)
      {
        if (param_.sigma_[k])
        { sum += Law::Normal::lpdf(p_data()->elt(i,j), param_.mean_[k][j], param_.sigma_[k]);}
      }
      return sum;
    }
    /** Initialize randomly the parameters of the Gaussian mixture. The centers
     *  will be selected randomly among the data set and the standard-deviations
     *  will be set to 1.
     */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** Compute the weighted mean and the common standard deviation. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return this->nbCluster()*this->nbVariable() + this->nbCluster();}
};

/* Initialize randomly the parameters of the Gaussian mixture. The centers
 *  will be selected randomly among the data set and the standard-deviation
 *  will be set to 1.
 */
template<class Array>
void MixtureGaussian_sk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  this->randomMean(p_tik);
  // compute the standard deviation
  Real variance;
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    variance = sqrt( ( p_tik->col(k).transpose()
                     *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
                      ).square()
                     ).sum() / (p_data()->sizeCols()*p_nk->elt(k))
                   );
    param_.sigma_[k] = ((variance<=0) || !Arithmetic<Real>::isFinite(variance))
                       ? 1.
                       : std::sqrt(variance/(this->nbSample()*this->nbVariable()));
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("MixtureGaussian_sk<Array>::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk)  done\n");
#endif
}

/* Compute the weighted mean and the common standard deviation. */
template<class Array>
bool MixtureGaussian_sk<Array>::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  // compute the means
  if (!this->updateMean(p_tik)) return false;
  // compute the standard deviation
  for (int k= p_tik->beginCols(); k < p_tik->endCols(); ++k)
  {
    param_.sigma_[k]
    = sqrt( ( p_tik->col(k).transpose()
             *(*p_data() - (Const::Vector<Real>(p_data()->rows()) * param_.mean_[k])
              ).square()
            ).sum()
           /(p_data()->sizeCols()*p_nk->elt(k))
          );
//    if (param(k).sigma_ <= 0.) return false;
  }
  return true;
}

} // namespace STK

#endif /* STK_MIXTUREGAUSSIAN_SK_H */
