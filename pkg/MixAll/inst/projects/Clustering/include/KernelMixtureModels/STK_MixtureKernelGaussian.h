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
 * created on: Oct 24, 2014
 * Author:   Serge Iovleff
 **/

/** @file STK_MixtureKernelGaussian.h
 *  @brief In this file we define the MixtureKernelGaussian_sk class
 **/

#ifndef STK_KERNELGAUSSIAN_H
#define STK_KERNELGAUSSIAN_H

#include <Arrays/include/STK_Array2DPoint.h>
#include <STatistiK/include/STK_Stat_Online.h>
#include "../STK_IMixtureDensity.h"

namespace STK
{

//forward declaration, to allow for recursive template
class MixtureKernelGaussian_sk;
class MixtureKernelGaussian_s;

namespace hidden
{
/** @ingroup Clustering
 *  Traits class for the MixtureKernelGaussian_sk traits policy. */
template<>
struct MixtureTraits< MixtureKernelGaussian_sk >
{
  typedef ArrayXX Array;
  typedef typename Array::Type       Type;
  /** Type of the structure storing the parameters of a MixtureKernelGaussian_sk model*/
  typedef ModelParameters<Clust::KernelGaussian_sk_> Parameters;
};

/** @ingroup Clustering
 *  Traits class for the MixtureKernelGaussian_sk traits policy. */
template<>
struct MixtureTraits< MixtureKernelGaussian_s >
{
  typedef ArrayXX Array;
  typedef typename Array::Type       Type;
  /** Type of the structure storing the parameters of a MixtureKernelGaussian_s model*/
  typedef ModelParameters<Clust::KernelGaussian_s_> Parameters;
};

} // namespace Clust

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixtureKernelGaussian_sk model.
 */
template<>
struct ModelParameters<Clust::KernelGaussian_sk_>
{
    /** variance of the variables */
    CPointX sigma2_;
    /** dimension of the gaussian kernel */
    CPointX dim_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): sigma2_(nbCluster, 1.), dim_(nbCluster, 2) {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : sigma2_(param.sigma2_), dim_(param.dim_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the standard deviation of the kth cluster */
    inline Real const& sigma2(int k) const { return sigma2_[k];}
    /** @return the dimension of the kth cluster */
    inline Real const& dim(int k) const { return dim_[k];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    {
      for (int k = sigma2_.begin(); k< sigma2_.end(); ++k)
      {
        sigma2_[k] = 1.;
        dim_[k] = 2.;
      }
    }
};

/** @ingroup Clustering
 *  Structure encapsulating the parameters of a MixtureKernelGaussian_s model.
 */
template<>
struct ModelParameters<Clust::KernelGaussian_s_>
{
    /** variance of the variables */
    Real sigma2_;
    /** dimension of the gaussian kernel */
    CPointX dim_;
    /** default constructor
     *  @param nbCluster the number of class of the mixture
     **/
    ModelParameters(int nbCluster): sigma2_(1.), dim_(nbCluster, 2.) {}
    /** copy constructor.
     *  @param param the parameters to copy.
     **/
    ModelParameters( ModelParameters const& param)
                   : sigma2_(param.sigma2_), dim_(param.dim_)
    {}
    /** destructor */
    ~ModelParameters() {}
    /** @return the standard deviation of the kth cluster */
    inline Real const& sigma2(int k) const { return sigma2_;}
    /** @return the dimension of the kth cluster */
    inline Real const& dim(int k) const { return dim_[k];}
    /** resize the set of parameter */
    inline void resize(Range const& range)
    { sigma2_ = 1.; dim_ = 2.;}
};

/** @ingroup Clustering
 *  Base class for the Kernel models
 **/
template<class Derived>
struct KernelHandlerBase: public IRecursiveTemplate<Derived>
{
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& sigma2(int k) const { return this->asDerived().sigma2Impl(k);}
  /** @return the value of lambda of the kth cluster and jth variable */
  inline Real const& dim(int k) const { return this->asDerived().dimImpl(k);}
};

/** @ingroup Clustering
 *  The Gaussian mixture model @c MixtureKernelGaussian_sk is an isotrope Gaussian
 *  mixture model on a kernel space. It has a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k
 *    \sum_{k=1}^K p_k \left(\frac{1}{\sqrt{2\pi}\sigma_k}\right)^{d_k}
 *    \exp\left\{ -\frac{\|\phi(x)-m_k\|^2}{2\sigma_k^2} \right\}
 * \f]
 * where \f$ \phi \f$ denote a feature mapping from the original space to an RKHS.
 *
 * In a MixtureKernelGaussian_sk model, the data set refer to the Gram's matrix.
 **/
class MixtureKernelGaussian_sk : public IMixtureDensity<MixtureKernelGaussian_sk >
{
  public:
    typedef IMixtureDensity<MixtureKernelGaussian_sk > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline MixtureKernelGaussian_sk( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixtureKernelGaussian_sk( MixtureKernelGaussian_sk const& model): Base(model) {}
    /** destructor */
    inline ~MixtureKernelGaussian_sk() {}
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const
    { return param_.dim_.sum() + this->nbCluster();}
    /** set the dimensions of the kernel mixture model using an unique value */
    inline void setDim(Real const& dim)  { param_.dim_ = dim;}
    /** set the dimension of the kernel mixture model */
    inline void setDim(PointX const& dim)  { param_.dim_ = dim;}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      return(- p_data()->elt(i,k)/(2.*param_.sigma2_[k])
             - (std::log(param_.sigma2_[k])+2.*Const::_LNSQRT2PI_)*param_.dim_[k]/2.);
    }
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param pk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    Real impute(int i, int j, Weights const& pk) const { return 0.;}
    /** @return a simulated value for the jth variable of the ith sample */
    inline Real rand(int i, int j, int k) const { return 0.;}
    /** Initialize randomly the variances of the Gaussian kernel mixture. */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** update the variances. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
};

/** @ingroup Clustering
 *  The Gaussian mixture model @c MixtureKernelGaussian_sk is an isotrope Gaussian
 *  mixture model on a kernel space. It has a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k
 *    \sum_{k=1}^K p_k \left(\frac{1}{\sqrt{2\pi}\sigma_k}\right)^{d_k}
 *    \exp\left\{ -\frac{\|\phi(x)-m_k\|^2}{2\sigma_k^2}  \right\}
 * \f]
 * where \f$ \phi \f$ denote a feature mapping from the original space to an RKHS.
 *
 * In a MixtureKernelGaussian_sk model, the data set refer to the Gram's matrix.
 **/
class MixtureKernelGaussian_s : public IMixtureDensity<MixtureKernelGaussian_s >
{
  public:
    typedef IMixtureDensity<MixtureKernelGaussian_s > Base;
    using Base::param_;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    inline MixtureKernelGaussian_s( int nbCluster): Base(nbCluster) {}
    /** copy constructor
     *  @param model The model to copy
     **/
    inline MixtureKernelGaussian_s( MixtureKernelGaussian_s const& model): Base(model) {}
    /** destructor */
    inline ~MixtureKernelGaussian_s() {}
    /** @return a constant reference on the paremeter handler structure.*/
    inline Parameters const& getParameters() const { return param_;}
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return param_.dim_.sum() + 1;}
    /** set the dimensions of the kernel mixture model using an unique value */
    inline void setDim(Real const& dim)  { param_.dim_ = dim;}
    /** set the dimensions of the kernel mixture model */
    inline void setDim(PointX const& dim)  { param_.dim_ = dim;}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      return(- p_data()->elt(i,k)/(2.*param_.sigma2_)
             - (std::log(param_.sigma2_)+2.*Const::_LNSQRT2PI_)*param_.dim_[k]/2.);
    }
    /** @return an imputation value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param tk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    inline Real impute(int i, int j, Weights const& tk) const { return 0.;}
    /** @return a simulated value for the jth variable of the ith sample */
    inline Real rand(int i, int j, int k) const { return 0.;}
    /** Initialize randomly the variances of the Gaussian kernel mixture. */
    void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
    /** update the variances. */
    bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
};

/* Initialize randomly the parameters of the Gaussian mixture. */
inline void MixtureKernelGaussian_sk::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering MixtureKernelGaussian_sk::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) \n");
#endif
  param_.sigma2_ = sum( p_data()->prod(*p_tik) )/ (*p_nk * param_.dim_)
                 + CPointX(p_tik->cols()).rand(Law::Normal(0, 0.05)).abs();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("MixtureKernelGaussian_sk::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk)  done\n");
#ifdef STK_MIXTURE_DEBUG
  stk_cout << param_.sigma2_ << "\n";
#endif
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
inline bool MixtureKernelGaussian_sk::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("Entering MixtureKernelGaussian_sk::run( CArrayXX const*  p_tik, CPointX const* p_nk) \n");
#endif
  param_.sigma2_ =  sum( p_data()->prod(*p_tik) )/ (*p_nk * param_.dim_);
  //if ((param_.sigma2_() <= 0.).any()) return false; // not work with Array1D
#ifdef STK_MIXTURE_DEBUG
  stk_cout << param_.sigma2_ << "\n";
#endif
  return true;
}

/* Initialize randomly the parameters of the Gaussian mixture. */
inline void MixtureKernelGaussian_s::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  // compute the standard deviation
  param_.sigma2_ = p_data()->prod(*p_tik).sum()/(this->nbSample() * param_.dim_.sum())
                 + std::abs(Law::generator.randGauss(0, 0.05));
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("MixtureKernelGaussian_s::randomInit( CArrayXX const*  p_tik, CPointX const* p_nk)  done\n");
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
inline bool MixtureKernelGaussian_s::run( CArrayXX const*  p_tik, CPointX const* p_nk) 
{
  param_.sigma2_ =  ( p_data()->prod( *p_tik ) ).sum()/(this->nbSample() * param_.dim_.sum());
  if (param_.sigma2_ <= 0.)  return false;
  return true;
}

} // namespace STK

#endif /* STK_KERNELGAUSSIAN_H */
