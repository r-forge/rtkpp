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

/** @file STK_PoissonParametersHandler.h
 *  @brief In this file we define the parameters handlers classes for Poisson
 *  mixture models
 **/

#ifndef STK_POISSONPARAMETERSHANDLER_H
#define STK_POISSONPARAMETERSHANDLER_H

#include <STatistiK/include/STK_Stat_Online.h>
#include "../PoissonModels/STK_PoissonParameters.h"


namespace STK
{

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Poisson_ljlk model
 **/
template <>
struct ParametersHandler<Clust::Poisson_ljlk_>
{
    typedef ModelParameters<Clust::Poisson_ljlk_> Parameters;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<Real, Real> > stat_lambdak_;
    /** Array of the lambdaj_ statistics */
    Stat::Online<CVectorX, Real>  stat_lambdaj_;
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
    /** Initialize the statistics. */
    void resize(Range const& range);
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Poisson_ljk model
 **/
template <>
struct ParametersHandler<Clust::Poisson_ljk_>
{
    typedef ModelParameters<Clust::Poisson_ljk_> Parameters;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_lambda_;
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
    /** Initialize the statistics. */
    void resize(Range const& range);
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for Poisson_lk model
 **/
template <>
struct ParametersHandler<Clust::Poisson_lk_>
{
    typedef ModelParameters<Clust::Poisson_lk_> Parameters;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<Real, Real> > stat_lambda_;
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
    /** Initialize the statistics. */
    void resize(Range const& range);
};


} // namespace STK

#endif /* STK_POISSONPARAMETERSHANDLER_H */
