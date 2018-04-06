/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  stkpp::
 * created on: Sep 8, 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_KernelParametersHandler.h
 *  @brief In this file we define the parameters handlers classes for kernel
 *  mixture models.
 **/


#ifndef STK_KERNELPARAMETERSHANDLER_H
#define STK_KERNELPARAMETERSHANDLER_H

#include <STatistiK/include/STK_Stat_Online.h>
#include "STK_KernelParameters.h"

namespace STK
{
/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Kmm_s model
 **/
template <>
struct ParametersHandler<Clust::Kmm_s_>
{
    typedef ModelParameters<Clust::Kmm_s_> Parameters;
    /** sigma2 statistics */
    Stat::Online<Real, Real> stat_sigma2_;
    /** Array of the dim statistics */
    Array1D< Stat::Online<Real, Real> > stat_dim_;
    /** default constructor. All lambdas are initialized to 1. */
    ParametersHandler( int nbCluster);
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    ParametersHandler( ParametersHandler const& param);
    /** destructor */
    ~ParametersHandler();
    /** copy operator */
    ParametersHandler& operator=( ParametersHandler const& other);
    /** update statistics of the parameters. */
    void updateStatistics(Parameters const& param);
    /** Set the computed statistics */
    void setStatistics(Parameters& param);
    /** Release the computed statistics */
    void releaseStatistics();
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Kmm_sk model
 **/
template <>
struct ParametersHandler<Clust::Kmm_sk_>
{
    typedef ModelParameters<Clust::Kmm_sk_> Parameters;
    /** Array of the sigma2 statistics */
    Array1D< Stat::Online<Real, Real> > stat_sigma2_;
    /** Array of the dim statistics */
    Array1D< Stat::Online<Real, Real> > stat_dim_;
    /** default constructor. All lambdas are initialized to 1. */
    ParametersHandler( int nbCluster);
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    ParametersHandler( ParametersHandler const& param);
    /** destructor */
    ~ParametersHandler();
    /** copy operator */
    ParametersHandler& operator=( ParametersHandler const& other);
    /** update statistics of the parameters. */
    void updateStatistics(Parameters const& param);
    /** Set the computed statistics */
    void setStatistics(Parameters& param);
    /** Release the computed statistics */
    void releaseStatistics();
};

} // namespace STK
#endif /* STK_KERNELPARAMETERSHANDLER_H_ */
