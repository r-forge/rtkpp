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

/** @file learnMixture.cpp
 *  @brief In this file we launch the computation for estimating a mixture model.
 **/


#include "RTKpp.h"
#include "../inst/projects/MixAll/LearnLauncher.h"
#include "../inst/projects/MixAll/RDataHandler.h"

/** @param model ClusterDiagModel S4 class
 *  @param models a vector of string with the model names to try
 */
RcppExport SEXP learnMixture( SEXP model, SEXP models
                            , SEXP algo, SEXP critName
                            , SEXP nbCore);

/** @param model ClusterMixedData S4 class
 *  @param algo estimation strategy S4 class
 *  @param critName name criteria string
 *  @param nbCore number of core to use
 */
RcppExport SEXP learnMixedData( SEXP model
                              , SEXP algo, SEXP critName
                              , SEXP nbCore);

/** @param model a ClusterDiagModel S4 class
 *  @param models a vector of string with the model names to try
 */
RcppExport SEXP learnKernelMixture( SEXP model, SEXP models
                                  , SEXP algo, SEXP critName
                                  , SEXP nbCore);

