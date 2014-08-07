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

using namespace Rcpp;
using namespace STK;

/**  @param tab2d R matrix
 */
RcppExport SEXP wrapper( SEXP tab2d )
{
  BEGIN_RCPP

  // wrap SEXP
  NumericMatrix RData(tab2d);

  // wrap Rcpp matrix
  RcppMatrix<double> data(RData);
  for (int i=data.beginRows(); i< data.endRows(); ++i)
  {
    for (int j= data.beginCols(); j < data.endCols(); ++j)
    {
      if (Arithmetic<Real>::isNA(data(i,j)))
      { data(i,j)= 100*i+j;}
    }
  }
  data(0,0) = Arithmetic<double>::NA();
  // return final output
  return tab2d;

  END_RCPP
}

