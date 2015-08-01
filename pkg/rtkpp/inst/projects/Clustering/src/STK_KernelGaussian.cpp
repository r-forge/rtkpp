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
 * created on: Oct 24, 2013
 * Author:   Serge Iovleff
 **/

/** @file STK_KernelGaussian.cpp
 *  @brief In this file we implement the KernelGaussian_s and KernelGaussian_sk classes
 **/

#include "../include/KernelMixtureModels/STK_KernelGaussian.h"
#include <Arrays/include/STK_Array2DVector.h>

namespace STK
{
/* */
void KernelGaussian_sk::getParameters(ArrayXX& param) const
{
  param.resize(2, nbCluster());
  for (int k= param_.sigma2_.begin(); k < param_.sigma2_.end(); ++k)
  {
    param(baseIdx  , k) = std::sqrt(param_.sigma2_[k]);
    param(baseIdx+1, k) = param_.dim_[k];
  }
}
/* Write the parameters on the output stream os */
void KernelGaussian_sk::writeParameters(ostream& os) const
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    os << _T("---> Component ") << k << _T("\n");
    os << _T("sigma = ") << std::sqrt(param_.sigma2_[k]) << _T("\n");
    os << _T("dim = ")   << param_.dim_[k] << _T("\n");
  }
}
/* Initialize randomly the parameters of the Gaussian mixture.
 */
void KernelGaussian_sk::randomInit()
{
  // compute the standard deviation
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  { param_.sigma2_[k] = p_data()->col(k).dot(p_tik()->col(k))/p_nk()->elt(k);}

#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelGaussian_sk::randomInit() done\n");
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
bool KernelGaussian_sk::mStep()
{
  param_.sigma2_ =  sum( p_data()->prod(*p_tik()) )/ (*p_nk() * param_.dim_);
  if ((param_.sigma2_ <= 0.).any())  return false;

  return true;
}

/* */
void KernelGaussian_s::getParameters(ArrayXX& param) const
{
  param.resize(2, nbCluster());
  for (int k= param_.dim_.begin(); k < param_.dim_.end(); ++k)
  {
    param(baseIdx  , k) = std::sqrt(param_.sigma2_);
    param(baseIdx+1, k) = param_.dim_[k];
  }
}
/* Write the parameters on the output stream os */
void KernelGaussian_s::writeParameters(ostream& os) const
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    os << _T("---> Component ") << k << _T("\n");
    os << _T("sigma = ") << std::sqrt(param_.sigma2_) << _T("\n");
    os << _T("dim = ")  << param_.dim_[k] << _T("\n");
  }
}
/* Initialize randomly the parameters of the Gaussian mixture.
 */
void KernelGaussian_s::randomInit()
{
  // compute the standard deviation
  param_.sigma2_ = 0.;
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  { param_.sigma2_ += p_data()->col(k).dot(p_tik()->col(k));}
  param_.sigma2_ /= this->nbSample();

#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelGaussian_s::randomInit() done\n");
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
bool KernelGaussian_s::mStep()
{
  param_.sigma2_ =  ( p_data()->prod( *p_tik() ) ).sum()/ (this->nbSample() * param_.dim_.sum());
  if (param_.sigma2_ <= 0.)  return false;

  return true;
}

} // namespace STK
