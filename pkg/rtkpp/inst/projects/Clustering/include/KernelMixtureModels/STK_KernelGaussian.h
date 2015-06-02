/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015 Serge Iovleff

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

/** @file STK_KernelGaussian.h
 *  @brief In this file we define the KernelGaussian_sk class
 **/

#ifndef STK_KERNELGAUSSIAN_H
#define STK_KERNELGAUSSIAN_H

#include <Arrays/include/STK_Array2DVector.h>
#include <STatistiK/include/STK_Stat_Online.h>
#include "../STK_IMixtureModel.h"

namespace STK
{

//forward declaration, to allow for recursive template
class KernelGaussian_sk;

namespace Clust
{
/** @ingroup Clustering
 *  Traits class for the KernelGaussian_sk traits policy. */
template<>
struct MixtureTraits< KernelGaussian_sk >
{
  typedef ArrayXX Array;
  typedef typename Array::Type       Type;
  typedef Array2D<Real>              Param;
  typedef ParametersHandler<KernelGaussian_sk_>  ParamHandler;
};

} // namespace Clust

/** Specialization of the ParametersHandler struct for KernelGaussian_sk models*/
template <>
struct ParametersHandler<Clust::KernelGaussian_sk_>
{
  /** default constructor */
  inline ParametersHandler(int nbCluster): sigma2_(nbCluster, 1.), stat_sigma2_(nbCluster)
                                         , lambda_(nbCluster, 1.), stat_lambda_(nbCluster)
  {}
  /** copy constructor.
   * @param param the parameters to copy.
   **/
  inline ParametersHandler( ParametersHandler const& param)
                          : sigma2_(param.sigma2_)
                          , stat_sigma2_(param.stat_sigma2_)
                          , lambda_(param.lambda_)
                          , stat_lambda_(param.stat_lambda_)
  {}
  /** destructor */
  inline ~ParametersHandler() {}
  /** Store the intermediate results of the Mixture.
   *  @param iteration Provides the iteration number beginning after the burn-in period.
   **/
  inline void storeIntermediateResults(int iteration)
  { stat_sigma2_.update(sigma2_); stat_lambda_.update(lambda_);}
  /** Release the stored results. This is usually used if the estimation
   *  process failed.
   **/
  inline void releaseIntermediateResults()
  { stat_sigma2_.release(); stat_lambda_.release();}
  /** set the parameters stored in stat_proba_ and release stat_proba_. */
  inline void setParameters()
  {
    sigma2_ = stat_sigma2_.mean_;
    stat_sigma2_.release();
    lambda_ = stat_lambda_.mean_;
    stat_lambda_.release();
  }
  /** vector of the standard deviations */
  VectorX sigma2_;
  /** Array of the statistics */
  Stat::Online<VectorX, Real> stat_sigma2_;
  /** vector of the dimensions */
  VectorX lambda_;
  /** Array of the statistics */
  Stat::Online<VectorX, Real> stat_lambda_;
};


/** @ingroup Clustering
 *  The Gaussian mixture model @c KernelGaussian_sk is
 *  an isotrope Gaussian model on a kernel space. It has a density function of the
 *  form
 * \f[
 *  f(\mathbf{x}|\theta) = \sum_{k=1}^K p_k \prod_{j=1}^d
 *    \frac{1}{\sqrt{2\pi}\sigma^j_{k}} \exp\left\{-\frac{(\phi(x^j_)-\mu^j_{ik})^2}{2(\sigma^j_{k})^2}\right\}.
 * \f]
 * The data set
 **/
class KernelGaussian_sk : public IMixtureModel<KernelGaussian_sk >
{
  public:
    typedef IMixtureModel<KernelGaussian_sk > Base;
    typedef ParametersHandler<Clust::KernelGaussian_sk_> ParamHandler;

    using Base::p_tik; using Base::param_;
    using Base::p_nk;
    using Base::p_data;

    /** default constructor
     * @param nbCluster number of cluster in the model
     **/
    KernelGaussian_sk( int nbCluster);
    /** copy constructor
     *  @param model The model to copy
     **/
    KernelGaussian_sk( KernelGaussian_sk const& model);
    /** destructor */
    ~KernelGaussian_sk();
    /** @return a constant reference on the paremeter handler structure.*/
    inline ParamHandler const& getParameters() const { return param_;}
    /** @return the number of free parameters of the model */
    inline int computeNbFreeParameters() const { return this->nbCluster();}
    /** set the dimension of the kernel mixture model */
    inline void setLambda(Real const& lambda)  { param_.lambda_ = lambda;}
    /** set the dimension of the kernel mixture model */
    inline void setLambda(VectorX const& lambda)  { param_.lambda_ = lambda;}
    /** @return the value of the probability of the i-th sample in the k-th component.
     *  @param i,k indexes of the sample and of the component
     **/
    inline Real lnComponentProbability(int i, int k) const
    {
      return(- p_data()->elt(i,k)/(2.*param_.sigma2_[k])
             - (std::log(param_.sigma2_[k])+2.*Const::_LNSQRT2PI_)*param_.lambda_[k]/2.);
    }
    /** get the parameters */
    void getParameters(ArrayXX& param) const;
    /** Write the parameters on the output stream os */
    void writeParameters(ostream& os) const;
    /** Initialize randomly the variances of the Gaussian kernel mixture. */
    void randomInit();
    /** update the variances. */
    bool mStep();
};

} // namespace STK

#endif /* STK_KERNELGAUSSIAN_H */
