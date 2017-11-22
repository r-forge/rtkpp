/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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
 * Project:  stkpp::Arrays
 * created on: 17 févr. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Arrays_Util.h
 *  @brief In this file we define utilities functions and enum for the Array classes.
 **/


#ifndef STK_ARRAY_UTIL_H
#define STK_ARRAY_UTIL_H

#include "STKernel/include/STK_Range.h"

namespace STK
{

namespace Arrays
{

/** @ingroup Arrays
 *  Intrinsic dimension of the container : 1D, 2D, 3D or 4D. 0D is for scalar
 **/
enum Dimension
{
  _0D_ = 0, ///< a single scalar have no dimension
  _1D_ = 1,
  _2D_ = 2,
  _3D_ = 3,
  _4D_ = 4
};
/** @ingroup Arrays
 *  Define the Storage Orientation of the container
 **/
enum Orientation
{
  by_row_ =0,  ///< storage by row
  by_col_ =1   ///< storage by column
};

/** @ingroup Arrays
 *  Define the different type of Array that can be handle by STK++
 **/
enum Storage
{
  dense_ =1,  ///< dense matrix/vector/array/expression
  sparse_=0   ///< sparse matrix/vector/array/expression
};

/**  @ingroup Arrays
 *   structures of Arrays that can be handled by STK++
 **/
enum Structure
{
  array2D_ =0 ,       ///< general matrix/array/expression
  square_,            ///< square matrix/array/expression
  diagonal_,          ///< diagonal matrix/array/expression
  lower_triangular_,  ///< lower triangular matrix/array/expression
  upper_triangular_,  ///< upper triangular matrix/array/expression
  lower_symmetric_,    ///< lower symmetric matrix/array/expression
  upper_symmetric_,    ///< upper symmetric matrix/array/expression
  vector_,            ///< column oriented vector/array/expression
  point_,             ///< row oriented vector/array/expression
  number_,            ///< (1,1) matrix/vector/array/expression (like a number)
  expression_         ///< An expression that will be evaluated further
};

/** @ingroup Arrays
 *  @return n+m, where n is the first number such that m < 2^n.
 *  @param m the size of the container
 **/
inline int evalSizeCapacity(int m)
{
  int n = 0;
  for (int k=1; k <= m; k <<= 1) {n++;}
  return(m+n);
}


/** @ingroup Arrays
 *  @return range of size n+m, where n is the first number such that m < 2^n.
 *  @param I the range of the container
 **/
inline Range evalRangeCapacity(STK::Range const& I)
{
  int n = 0;
  for (int k=1; k <= I.size(); k <<= 1){ n++;}
  return Range(I.begin(),I.size() + n);
}

/** @ingroup Arrays
 *  convert an Arrays::Structure to a String.
 *  @param type the type of Structure we want to convert
 *  @return the string associated to this type.
 **/
inline std::string structureToString( Structure const& type)
{
  if (type == array2D_)          return String(_T("array2D_"));
  if (type == square_)           return String(_T("square_"));
  if (type == diagonal_)         return String(_T("diagonal_"));
  if (type == lower_triangular_) return String(_T("lower_triangular_"));
  if (type == upper_triangular_) return String(_T("upper_triangular_"));
  if (type == vector_)           return String(_T("vector_"));
  if (type == point_)            return String(_T("point_"));
  if (type == number_)           return String(_T("number_"));
  if (type == expression_)       return String(_T("expression_"));
  return String(_T("unknown"));
}


} // namespace Arrays

} // namespace STK



#endif /* STK_ARRAY_UTIL_H */
