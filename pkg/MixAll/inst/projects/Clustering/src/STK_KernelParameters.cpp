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

/** @file STK_KernelParameters.cpp
 *  @brief In this file we implement the ModelParameters class for kernel
 *  mixture models
 **/


#include <Clustering/include/KernelModels/STK_KernelParameters.h>

namespace STK
{
/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Kmm_s_>::ModelParameters(int nbCluster)
                             : sigma2_(1.), dim_(nbCluster, 10) {}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Kmm_s_>::ModelParameters( ModelParameters const& param)
               : sigma2_(param.sigma2_), dim_(param.dim_)
{}
/* destructor */
ModelParameters<Clust::Kmm_s_>::~ModelParameters() {}



/* default constructor
 *  @param nbCluster the number of class of the mixture
 **/
ModelParameters<Clust::Kmm_sk_>::ModelParameters(int nbCluster)
                             : sigma2_(nbCluster, 1.), dim_(nbCluster, 10) {}
/* copy constructor.
 *  @param param the parameters to copy.
 **/
ModelParameters<Clust::Kmm_sk_>::ModelParameters( ModelParameters const& param)
               : sigma2_(param.sigma2_), dim_(param.dim_)
{}
/* destructor */
ModelParameters<Clust::Kmm_sk_>::~ModelParameters() {}


} // namespace STK
