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
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_GammaParametersHandler.h
 *  @brief In this file we define the parameters handlers classes for gamma
 *  mixture models
 **/

#ifndef STK_GAMMAPARAMETERSHANDLER_H
#define STK_GAMMAPARAMETERSHANDLER_H

#include <STatistiK/include/STK_Stat_Online.h>
#include "../GammaModels/STK_GammaParameters.h"

namespace STK
{
/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_a_bjk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_a_bjk_>
{
    typedef ModelParameters<Clust::Gamma_a_bjk_> Parameters;
    /** shape statistics */
    Stat::Online<Real, Real> stat_shape_;
    /** Array of the scale statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_a_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_a_bk_>
{
    typedef ModelParameters<Clust::Gamma_a_bk_> Parameters;
    /** shape statistics */
    Stat::Online<Real, Real> stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
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
    /** Set the computed statistics */
    void releaseStatistics();
    /** Initialize the statistics. */
    void resize(Range const& range);
};


/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_aj_bjk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_aj_bjk_>
{
    typedef ModelParameters<Clust::Gamma_aj_bjk_> Parameters;
    /** shape statistics */
    Stat::Online<CPointX, Real> stat_shape_;
    /** Array of the scale statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_aj_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_aj_bk_>
{
    typedef ModelParameters<Clust::Gamma_aj_bk_> Parameters;
    /** shape statistics */
    Stat::Online<CPointX, Real> stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_ajk_b model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_b_>
{
    typedef ModelParameters<Clust::Gamma_ajk_b_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_ajk_bj model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_bj_>
{
    typedef ModelParameters<Clust::Gamma_ajk_bj_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<CPointX, Real> stat_scale_;
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
 *  Specialization of the ParametersHandler struct for STK::Gamma_ajk_bjk
 *  mixture model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_bjk_>
{
    typedef ModelParameters<Clust::Gamma_ajk_bjk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_ajk_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ajk_bk_>
{
    typedef ModelParameters<Clust::Gamma_ajk_bk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_ak_b model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_b_>
{
    typedef ModelParameters<Clust::Gamma_ak_b_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<Real, Real> stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_ak_bj model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_bj_>
{
    typedef ModelParameters<Clust::Gamma_ak_bj_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Stat::Online<CPointX, Real> stat_scale_;
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
 *  Specialization of the ParametersHandler struct for Gamma_ak_bjk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_bjk_>
{
    typedef ModelParameters<Clust::Gamma_ak_bjk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_scale_;
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
    /** Set the computed statistics */
    void releaseStatistics();
    /** Initialize the statistics. */
    void resize(Range const& range);
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Gamma_ak_bk model
 **/
template <>
struct ParametersHandler<Clust::Gamma_ak_bk_>
{
    typedef ModelParameters<Clust::Gamma_ak_bk_> Parameters;
    /** Array of the mean statistics */
    Array1D< Stat::Online<Real, Real> > stat_shape_;
    /** Array of the standard deviation statistics */
    Array1D< Stat::Online<Real, Real> > stat_scale_;
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

#endif /* STK_GAMMAPARAMETERSHANDLER_H */
