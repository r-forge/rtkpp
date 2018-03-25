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

/** @file STK_CategoricalParametersHandler.h
 *  @brief In this file we define the parameters handlers classes for categorical
 *  mixture models
 **/

#ifndef STK_CATEGORICALPARAMETERSHANDLER_H
#define STK_CATEGORICALPARAMETERSHANDLER_H

#include <STatistiK/include/STK_Stat_Online.h>
#include "../CategoricalModels/STK_CategoricalParameters.h"

namespace STK
{

/** @ingroup Clustering
 * Specialization of the ParametersHandler struct for Categorical_pjk model */
template <>
struct ParametersHandler<Clust::Categorical_pjk_>
{
    typedef ModelParameters<Clust::Categorical_pjk_> Parameters;
    /** statistics of the probabilities */
    Array1D<  Stat::Online<CArrayXX, Real>  > stat_proba_;
    /** default constructor */
    ParametersHandler(int nbCluster);
    /** copy constructor */
    ParametersHandler(ParametersHandler const& model);
    /** copy operator */
    ParametersHandler& operator=( ParametersHandler const& other);
    /** update statistics of the parameters. */
    void updateStatistics(Parameters const& param);
    /** set and release the computed statistics */
    void setStatistics(Parameters& param);
    /** Set the computed statistics */
    void releaseStatistics();
    /** Initialize the parameters of the model. */
    void resize(Range const& rangeModalities, Range const& rangeData );
};

/** Specialization of the ParametersHandler struct for Categorical_pk model */
template <>
struct ParametersHandler<Clust::Categorical_pk_>
{
    typedef ModelParameters<Clust::Categorical_pk_> Parameters;
    /** statistics of the probabilities */
    Array1D<  Stat::Online<CVectorX, Real>  > stat_proba_;
    /** default constructor */
    ParametersHandler(int nbCluster);
    /** copy constructor */
    ParametersHandler(ParametersHandler const& model);
    /** copy operator */
    ParametersHandler& operator=( ParametersHandler const& other);
    /** update statistics of the parameters. */
    void updateStatistics(Parameters const& param);
    /** set the computed statistics */
    void setStatistics(Parameters& param);
    /** release the computed statistics */
    void releaseStatistics();
    /** Initialize the parameters of the model. */
    void resize(Range const& rangeModalities, Range const& rangeData );
};


} // namespace STK

#endif /* STK_CATEGORICALPARAMETERSHANDLER_H */
