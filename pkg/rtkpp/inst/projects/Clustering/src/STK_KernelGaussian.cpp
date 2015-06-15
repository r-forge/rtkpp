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
 *  @brief In this file we implement the KernelGaussian_sk class
 **/

#include "../include/KernelMixtureModels/STK_KernelGaussian.h"

namespace STK
{
/* default constructor
 * @param nbCluster number of cluster in the model
 **/
KernelGaussian_sk::KernelGaussian_sk( int nbCluster) : Base(nbCluster)
{}
/** copy constructor
 *  @param model The model to copy
 **/
KernelGaussian_sk::KernelGaussian_sk( KernelGaussian_sk const& model) : Base(model) {}
/** destructor */
KernelGaussian_sk::~KernelGaussian_sk() {}
/** */
void KernelGaussian_sk::getParameters(ArrayXX& param) const
{
  param.resize(nbCluster(), 2);
  for (int k= 0; k < nbCluster(); ++k)
  {
    param(baseIdx+k, baseIdx)   = std::sqrt(param_.sigma2_[baseIdx+k]);
    param(baseIdx+k, baseIdx+1) = param_.lambda_[baseIdx+k];
  }
}
/* Write the parameters on the output stream os */
void KernelGaussian_sk::writeParameters(ostream& os) const
{
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    os << _T("---> Component ") << k << _T("\n");
    os << _T("sigma2 = ") << param_.sigma2_[k] << _T("\n");
    os << _T("lambda = ") << param_.lambda_[k] << _T("\n");
  }
}
/* Initialize randomly the parameters of the Gaussian mixture.
 */
void KernelGaussian_sk::randomInit()
{
  // compute the standard deviation
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.sigma2_[k] = p_data()->col(k).dot(p_tik()->col(k))/p_nk()->elt(k);
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelGaussian<Array>::randomInit() done\n");
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
bool KernelGaussian_sk::mStep()
{
  // compute the standard deviation
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    param_.sigma2_[k] =  p_data()->col(k).dot(p_tik()->col(k))
                      / (p_nk()->elt(k)*param_.lambda_[k]);
    if (param_.sigma2_[k] <= 0.) return false;
  }
  return true;
}

} // namespace STK
