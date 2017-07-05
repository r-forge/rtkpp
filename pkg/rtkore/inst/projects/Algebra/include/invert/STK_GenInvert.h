/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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
 * Project:  stkpp::Algebra
 * created on: 9 août 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_GenInvert.h
 *  @brief In this file we implement inversion method for fixed size matrices.
 **/

#ifndef STK_GENINVERT_H
#define STK_GENINVERT_H

#ifdef STKUSELAPACK
#include <Algebra/include/STK_lapack_SymEigen.h>
#else
#include <Algebra/include/STK_SymEigen.h>
#endif

namespace STK
{

namespace hidden
{
/** @ingroup Algebra
 *  @brief compute the inverse of the 1x1 matrix and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertMatrix11( Matrix const& m, CArraySquare<typename Matrix::Type, 1>& inv)
{
  typedef typename Matrix::Type Type;
  // cofactor (0,0) [0]
  inv(0, 0) = Type(1);
  // compute determinant
  Type det = m(0,0);
  if (det == Type(0)) return Type(0);
  // compute inverse matrix
  inv /= det;
  return det;
}
/** @ingroup Algebra
 *  @brief compute the inverse of the 2x2 matrix m and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertMatrix22( Matrix const& m, CArraySquare<typename Matrix::Type, 2>& inv)
{
  typedef typename Matrix::Type Type;
  // cofactor (0,0) [0]
  inv(0, 0) =   m(1,1);
  // cofactor (1,0) [1]
  inv(1, 0) = - m(1, 0);
  // cofactor (1,1) [2]
  inv(0, 1) = - m(0,1);
  // cofactor (1,1) [3]
  inv(1, 1) =   m(0,0);

  // determinant
  Type det = m(0,0) * inv(0,0) + m(0,1) * inv(1,0);
  if (det == Type(0)) return Type(0);
  // compute inverse matrix
  inv /= det;
  return det;
}
/** @ingroup Algebra
 *  @brief compute the inverse of the 3x3 matrix m and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertMatrix33( Matrix const& m, CArraySquare<typename Matrix::Type, 3>& inv)
{
  typedef typename Matrix::Type Type;

  // cofactor
  inv(0,0) = m(1,1) * m(2,2) - m(2,1) * m(1,2);
  inv(0,1) = m(0,2) * m(2,1) - m(0,1) * m(2,2);
  inv(0,2) = m(0,1) * m(1,2) - m(0,2) * m(1,1);
  inv(1,0) = m(1,2) * m(2,0) - m(1,0) * m(2,2);
  inv(1,1) = m(0,0) * m(2,2) - m(0,2) * m(2,0);
  inv(1,2) = m(1,0) * m(0,2) - m(0,0) * m(1,2);
  inv(2,0) = m(1,0) * m(2,1) - m(2,0) * m(1,1);
  inv(2,1) = m(2,0) * m(0,1) - m(0,0) * m(2,1);
  inv(2,2) = m(0,0) * m(1,1) - m(1,0) * m(0,1);

  // computes the inverse of a matrix m
  Type det = m(0, 0) * inv(0,0) + m(0, 1) * inv(1,0) + m(0, 2) * inv(2,0);
  if (det == Type(0)) return Type(0);
  // compute inverse matrix
  inv /= det;
  return det;
}
/** @ingroup Algebra
 *  @brief compute the inverse of the 4x4 matrix m and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertMatrix44( Matrix const& m, CArraySquare<typename Matrix::Type, 4>& inv)
{
  typedef typename Matrix::Type Type;
  // cofactor (0,0) [0]
  inv(0,0) = m(1,1) * m(2,2) * m(3,3) - m(1,1) * m(3,2) * m(2,3)
           - m(1,2) * m(2,1) * m(3,3) + m(1,2) * m(3,1) * m(2,3)
           + m(1,3) * m(2,1) * m(3,2) - m(1,3) * m(3,1) * m(2,2);
  // cofactor (0,1) [4]
  inv(0,1) = -m(0,1) * m(2,2) * m(3,3) + m(0,1) * m(3,2) * m(2,3)
           +  m(0,2) * m(2,1) * m(3,3) - m(0,2) * m(3,1) * m(2,3)
           -  m(0,3) * m(2,1) * m(3,2) + m(0,3) * m(3,1) * m(2,2);
  // cofactor (0,2) [8]
  inv(0,2) = m(0,1) * m(1,2) * m(3,3) - m(0,1) * m(3,2) * m(1,3)
           - m(0,2) * m(1,1) * m(3,3) + m(0,2) * m(3,1) * m(1,3)
           + m(0,3) * m(1,1) * m(3,2) - m(0,3) * m(3,1) * m(1,2);
  // cofactor (0,3) [12]
  inv(0,3) = -m(0,1) * m(1,2) * m(2,3) + m(0,1) * m(2,2) * m(1,3)
           +  m(0,2) * m(1,1) * m(2,3) - m(0,2) * m(2,1) * m(1,3)
           -  m(0,3) * m(1,1) * m(2,2) + m(0,3) * m(2,1) * m(1,2);
  // cofactor (1,0) [1]
  inv(1,0) = -m(1,0) * m(2,2) * m(3,3) + m(1,0) * m(3,2) * m(2,3)
           +  m(1,2) * m(2,0) * m(3,3) - m(1,2) * m(3,0) * m(2,3)
           -  m(1,3) * m(2,0) * m(3,2) + m(1,3) * m(3,0) * m(2,2);
  // cofactor (1,1) [5]
  inv(1,1) = m(0,0) * m(2,2) * m(3,3) - m(0,0) * m(3,2) * m(2,3)
           - m(0,2) * m(2,0) * m(3,3) + m(0,2) * m(3,0) * m(2,3)
           + m(0,3) * m(2,0) * m(3,2) - m(0,3) * m(3,0) * m(2,2);
  // cofactor (1,2) [9]
  inv(1,2) = -m(0,0) * m(1,2) * m(3,3) + m(0,0) * m(3,2) * m(1,3)
           +  m(0,2) * m(1,0) * m(3,3) - m(0,2) * m(3,0) * m(1,3)
           -  m(0,3) * m(1,0) * m(3,2) + m(0,3) * m(3,0) * m(1,2);
  // cofactor (1,3) [13]
  inv(1,3) = m(0,0) * m(1,2) * m(2,3) - m(0,0) * m(2,2) * m(1,3)
           - m(0,2) * m(1,0) * m(2,3) + m(0,2) * m(2,0) * m(1,3)
           + m(0,3) * m(1,0) * m(2,2) - m(0,3) * m(2,0) * m(1,2);
  // cofactor (2,0) [2]
  inv(2,0) = m(1,0) * m(2,1) * m(3,3) - m(1,0) * m(3,1) * m(2,3)
           - m(1,1) * m(2,0) * m(3,3) + m(1,1) * m(3,0) * m(2,3)
           + m(1,3) * m(2,0) * m(3,1) - m(1,3) * m(3,0) * m(2,1);
  // cofactor (2,1)
  inv(2,1) = -m(0,0) * m(2,1) * m(3,3) + m(0,0) * m(3,1) * m(2,3)
           +  m(0,1) * m(2,0) * m(3,3) - m(0,1) * m(3,0) * m(2,3)
           -  m(0,3) * m(2,0) * m(3,1) + m(0,3) * m(3,0) * m(2,1);
  // cofactor (2,2)
  inv(2,2) = m(0,0) * m(1,1) * m(3,3) - m(0,0) * m(3,1) * m(1,3)
           - m(0,1) * m(1,0) * m(3,3) + m(0,1) * m(3,0) * m(1,3)
           + m(0,3) * m(1,0) * m(3,1) - m(0,3) * m(3,0) * m(1,1);
  // cofactor (2,3)
  inv(2,3) = -m(0,0) * m(1,1) * m(2,3) + m(0,0) * m(2,1) * m(1,3)
           +  m(0,1) * m(1,0) * m(2,3) - m(0,1) * m(2,0) * m(1,3)
           -  m(0,3) * m(1,0) * m(2,1) + m(0,3) * m(2,0) * m(1,1);
  // cofactor (3,0)
  inv(3,0) = -m(1,0) * m(2,1) * m(3,2) + m(1,0) * m(3,1) * m(2,2)
           +  m(1,1) * m(2,0) * m(3,2) - m(1,1) * m(3,0) * m(2,2)
           -  m(1,2) * m(2,0) * m(3,1) + m(1,2) * m(3,0) * m(2,1);
  // cofactor (3,1) [7]
  inv(3,1) = m(0,0) * m(2,1) * m(3,2) - m(0,0) * m(3,1) * m(2,2)
           - m(0,1) * m(2,0) * m(3,2) + m(0,1) * m(3,0) * m(2,2)
           + m(0,2) * m(2,0) * m(3,1) - m(0,2) * m(3,0) * m(2,1);
  // cofactor (3,2) [11]
  inv(3,2) = -m(0,0) * m(1,1) * m(3,2) + m(0,0) * m(3,1) * m(1,2)
           +  m(0,1) * m(1,0) * m(3,2) - m(0,1) * m(3,0) * m(1,2)
           -  m(0,2) * m(1,0) * m(3,1) + m(0,2) * m(3,0) * m(1,1);
  // cofactor (3,3) [15]
  inv(3,3) = m(0,0) * m(1,1) * m(2,2) - m(0,0) * m(2,1) * m(1,2)
           - m(0,1) * m(1,0) * m(2,2) + m(0,1) * m(2,0) * m(1,2)
           + m(0,2) * m(1,0) * m(2,1) - m(0,2) * m(2,0) * m(1,1);
  // compute determinant and inverse matrix
  Type det = m(0,0) * inv(0,0) + m(1,0) * inv(0,1) + m(2,0) * inv(0,2) + m(3,0) * inv(0,3);
  if (det == Type(0)) return Type(0);
  // compute inverse matrix
  inv /= det;
  return det;
}
/** @ingroup hidden utility class allowing to call the correct static function
 *  computing the inverse of a matrix. */
template<class Matrix, int Size> struct invertMatrixDispatcher;

template<class Matrix>
struct invertMatrixDispatcher<Matrix, 1>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 1>& inv)
  { return invertMatrix11(m, inv);}
};

template<class Matrix>
struct invertMatrixDispatcher<Matrix, 2>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 2>& inv)
  { return invertMatrix22(m, inv);}
};

template<class Matrix>
struct invertMatrixDispatcher<Matrix, 3>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 3>& inv)
  { return invertMatrix33(m, inv);}
};

template<class Matrix>
struct invertMatrixDispatcher<Matrix, 4>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 4>& inv)
  { return invertMatrix44(m, inv);}
};

} // namespace hidden

} // namespace STK

#endif /* STK_GENINVERT_H */
