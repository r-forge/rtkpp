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

/** @file STK_KernelGaussian_sk.cpp
 *  @brief In this file we implement the KernelGaussian_sk class
 **/


#include "../include/KernelMixtureModels/STK_KernelGaussian_sk.h"
#include <Analysis/include/STK_Const_Math.h>
#include <STatistiK/include/STK_Law_Normal.h>
#include <STatistiK/include/STK_Stat_Functors.h>

#if STK_MIXTURE_DEBUG | STK_MIXTURE_VERBOSE | STK_MIXTURE_VERY_VERBOSE
#include <Arrays/include/STK_Display.h>
#endif

namespace STK
{
/* default constructor
 * @param nbCluster number of cluster in the model
 **/
KernelGaussian_sk::KernelGaussian_sk( int nbCluster): Base(nbCluster) {}
/* copy constructor
 *  @param model The model to copy
 **/
KernelGaussian_sk::KernelGaussian_sk( KernelGaussian_sk const& model): Base(model) {}
/* destructor */
KernelGaussian_sk::~KernelGaussian_sk() {}
/* @return the number of free parameters of the model */
int KernelGaussian_sk::computeNbFreeParameters() const
{ return param_.dim_.sum() + this->nbCluster();}

/* @return the value of the probability of the i-th sample in the k-th component.
 *  @param i,k indexes of the sample and of the component
 **/
Real KernelGaussian_sk::lnComponentProbability(int i, int k) const
{
  return(- dik_.elt(i,k)/(2.*param_.sigma2_[k])
         - (std::log(param_.sigma2_[k])+2.*Const::_LNSQRT2PI_)*param_.dim_[k]/2.);
}

/* Initialize randomly the parameters of the Gaussian mixture. */
void KernelGaussian_sk::randomInit( CArrayXX const*  p_tik, CPointX const* p_tk)
{
#if STK_Kernel_DEBUG | STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering KernelGaussian_sk::randomInit( CArrayXX const*  p_tik, CPointX const* p_tk)\n");
#endif
  compute_dik(p_tik, p_tk);
  param_.sigma2_ = sum( dik_.prod(*p_tik) )/ (*p_tk * param_.dim_)
                 + CPointX(p_tik->cols()).rand(Law::Normal(0, 0.05)).abs();
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelGaussian_sk::randomInit( CArrayXX const*  p_tik, CPointX const* p_tk) done\n");
  stk_cout << param_.sigma2_ << "\n";
#endif
}

/* Compute the weighted means and the weighted standard deviations. */
bool KernelGaussian_sk::run( CArrayXX const*  p_tik, CPointX const* p_tk)
{
#if STK_Kernel_DEBUG | STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering KernelGaussian_sk::run( CArrayXX const*  p_tik, CPointX const* p_tk)\n");
#endif
  compute_dik(p_tik, p_tk);
#ifdef STK_Kernel_DEBUG
  stk_cout<< _T("Stat::sumByCol( dik_.prod(*p_tik) ) =\n")  << Stat::sumByCol( dik_.prod(*p_tik) ) << "\n";
#endif
  param_.sigma2_ =  Stat::sumByCol( dik_.prod(*p_tik) )/ (*p_tk * param_.dim_);
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("KernelGaussian_sk::run( CArrayXX const*  p_tik, CPointX const* p_tk) done\n");
  stk_cout << _T("sigma2 = ") << param_.sigma2_ << "\n";
#endif
  return true;
}

} // namespace STK
