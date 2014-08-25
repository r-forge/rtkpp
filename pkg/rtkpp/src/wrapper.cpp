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


#include "RTKpp.h"

using namespace Rcpp;
using namespace STK;

/**  @param tab2d R matrix
 */
RcppExport SEXP wrapper( SEXP tab1,SEXP tab2,SEXP tab3 )
{
  BEGIN_RCPP

  // wrap SEXP
  NumericMatrix RData1(tab1);
  NumericMatrix RData2(tab2);
  NumericMatrix RData3(tab3);

  // wrap Rcpp matrix
  RcppMatrix<double> data1(RData1);
  // start tests
  for (int i=data1.beginRows(); i< data1.endRows(); ++i)
  {
    for (int j= data1.beginCols(); j < data1.endCols(); ++j)
    {
      if (Arithmetic<Real>::isNA(data1(i,j)))
      { data1(i,j)= 100*i+j;}
    }
  }
  data1(0,0) = Arithmetic<double>::NA();
  RDataHandler handler(RData2, "model1", "gaussian_sjk");
  handler.addData(RData3,  "model2", "gaussian_sj");

  RcppMatrix<double> data2;
  int nbVar;
  handler.getData("model1", data2, nbVar);

  data2(0,0) = 10;

  Rcpp::List l;
  l.push_back(tab1, "1");
  data1(1,1) = 11;
  l.push_back(tab1, "2");
  data1(2,2) = 22;
  l.push_back(tab2, "3");

  // return final output
  return l;

  END_RCPP
}

