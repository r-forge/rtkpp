/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff

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
 * Project:  rtkpp::
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file wrap.cpp
 *  @brief In this file .
 **/


#include <Rcpp.h>
#include "../inst/include/rtkpp.h"
#include "../inst/include/DManager.h"

using namespace Rcpp;
using namespace STK;

/**  @param model ClusterDiagModel S4 class
 */
RcppExport SEXP clusterDiagGaussianModel( SEXP model, SEXP k, SEXP modelNames, SEXP strategy )
{
  BEGIN_RCPP

  // wrap S4 model with Rcpp
  S4 Rcppmodel(model);
  // wrap data matrix with Rcpp
  NumericMatrix RData = Rcppmodel.slot("data");
  // wrap Rcpp matrix with stk++
  RcppMatrix<double> data(RData);

  // start
  DataHandler handler;
  // ... TODO

  // return true
  return Rcpp::wrap(1);

  END_RCPP
}

