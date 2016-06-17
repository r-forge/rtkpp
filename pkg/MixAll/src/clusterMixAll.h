/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 28 juil. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file clusterMixture.cpp
 *  @brief In this file we launch the computation for estimating a mixture model.
 **/


#include "RTKpp.h"

/** @param model ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param models a vector of string with the model names to try
 */
RcppExport SEXP clusterMixture( SEXP model, SEXP nbCluster, SEXP models
                              , SEXP strategy, SEXP critName
                              , SEXP nbCore);

/** @param model ClusterMixedData S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param strategy estimation strategy S4 class
 *  @param critName name criteria string
 *  @param nbCore number of core to use
 */
RcppExport SEXP clusterMixedData( SEXP model, SEXP nbCluster
                                , SEXP strategy, SEXP critName
                                , SEXP nbCore);

/** @param model a ClusterDiagModel S4 class
 *  @param nbCluster a vector with the number of clusters to test
 *  @param models a vector of string with the model names to try
 */
RcppExport SEXP clusterKernelMixture( SEXP model, SEXP nbCluster, SEXP models
                                    , SEXP strategy, SEXP critName
                                    , SEXP nbCore);

/** Compute the Gram matrix and overwrite the data with the result.
 *  @param component a ClusterKernelComponent S4 class
 *  @param kernelName a string with the name of the kernel to use
 *  @param kernelParameters a vector with the optional parameters
 */
RcppExport SEXP clusterKernelCompute( SEXP component
                                    , SEXP kernelName
                                    , SEXP kernelParameters
                                    );

