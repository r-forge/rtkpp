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
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_DiagGaussianParametersHandler.h
 *  @brief In this file we define the bridge classes between the diagonal
 *  Gaussian mixtures and the composer.
 **/

#ifndef STK_DIAGGAUSSIANPARAMETERSHANDLER_H
#define STK_DIAGGAUSSIANPARAMETERSHANDLER_H

#include <STatistiK/include/STK_Stat_Online.h>
#include "../DiagGaussianModels/STK_DiagGaussianParameters.h"

namespace STK
{
/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for DiagGaussian_sjk model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_sjk_>
{
    typedef ModelParameters<Clust::Gaussian_sjk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_sigma_;
    /** default constructor. */
    ParametersHandler( int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
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
    /** Initialize the statistics. */
    void resize(Range const& range);
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for DiagGaussian_sk model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_sk_>
{
    typedef ModelParameters<Clust::Gaussian_sk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_sigma_;

    /** default constructor. */
    ParametersHandler( int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
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
    /** Initialize the statistics. */
    void resize(Range const& range);
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for DiagGaussian_sj model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_sj_>
{
    typedef ModelParameters<Clust::Gaussian_sj_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
   Stat::Online<CPointX, Real> stat_sigma_;

    /** default constructor. */
    ParametersHandler( int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
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
    /** Initialize the statistics. */
    void resize(Range const& range);
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for DiagGaussian_s model
 **/
template <>
struct ParametersHandler<Clust::Gaussian_s_>
{
    typedef ModelParameters<Clust::Gaussian_s_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_mean_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_sigma_;

    /** default constructor. */
    ParametersHandler( int nbCluster);
    /** copy constructor.
     *  @param param the parameters to copy.
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
    /** Initialize the statistics. */
    void resize(Range const& range);
};


} // namespace STK

#endif /* STK_DIAGGAUSSIANPARAMETERSHANDLER_H */
