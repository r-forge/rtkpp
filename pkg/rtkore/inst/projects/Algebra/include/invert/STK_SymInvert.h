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

/** @file STK_SymInvert.h
 *  @brief In this file we implement inversion method for fixed size symmetric matrices.
 **/

#ifndef STK_SYMINVERT_H
#define STK_SYMINVERT_H

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
 *  @brief compute the inverse of the symmetric matrix m of size 1x1 and store
 *  the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return The determinant value of m
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix11( Matrix const& m, CArraySquare<typename Matrix::Type, 1>& inv)
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
 *  @brief compute the inverse of the symmetric 2x2 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return The determinant value of m
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix22( Matrix const& m, CArraySquare<typename Matrix::Type, 2>& inv)
{
  typedef typename Matrix::Type Type;
  // cofactor (0,0) [0]
  inv(0, 0) =   m(1,1);
  // cofactor (1,0) [1]
  inv(1, 0) = - m(1, 0);
  // cofactor (1,1) [3]
  inv(1, 1) =   m(0,0);

  // symmetry
  inv(0,1) = inv(1,0);

  // determinant
  Type det = m(0,0) * inv(0,0) + m(1,0) * inv(1,0);
  if (det == Type(0)) return Type(0);
  // inverse matrix
  inv /= det;
  return det;
}
/** @ingroup Algebra
 *  @brief compute the inverse of the symmetric 3x3 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return The determinant value of m
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix33( Matrix const& m, CArraySquare<typename Matrix::Type, 3>& inv)
{
  typedef typename Matrix::Type Type;
  // cofactor (0,0) [0]
  inv(0, 0) = m(1,1) * m(2,2) - m(2, 1) * m(2, 1);
  // cofactor (1,0) [1]
  inv(1, 0) = m(1,2) * m(2,0) - m(1, 0) * m(2, 2);
  inv(2, 0) = m(1,0) * m(2,1) - m(2, 0) * m(1, 1);
  inv(1, 1) = m(0,0) * m(2,2) - m(0, 2) * m(2, 0);
  inv(2, 1) = m(2,0) * m(1,0) - m(0, 0) * m(2, 1);
  inv(2, 2) = m(0,0) * m(1,1) - m(1, 0) * m(0, 1);

  // symetric
  inv(0,1) = inv(1,0); inv(0,2) = inv(2,0);
  inv(1,2) = inv(2,1);

  // compute determinant and inverse matrix
  Type det = m(0,0) * inv(0,0) + m(1,0) * inv(1,0) + m(2,0) * inv(2,0);
  if (det == Type(0)) return Type(0);
  // inverse matrix
  inv /= det;
  return det;
}
/** @ingroup Algebra
 *  @brief compute the inverse of the symmetric 4x4 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return The determinant value of m
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix44( Matrix const& m, CArraySquare<typename Matrix::Type, 4>& inv)
{
  typedef typename Matrix::Type Type;
  // cofactor (0,0) [0]
  inv(0,0) = m(1,1) * m(2,2) * m(3,3) - m(1,1) * m(3,2) * m(3,2)
           - m(2,1) * m(2,1) * m(3,3) + m(2,1) * m(3,1) * m(3,2)
           + m(3,1) * m(2,1) * m(3,2) - m(3,1) * m(3,1) * m(2,2);
  // cofactor (1,0) [1]
  inv(1,0) = -m(1,0) * m(2,2) * m(3,3) + m(1,0) * m(3,2) * m(3,2)
           +  m(2,1) * m(2,0) * m(3,3) - m(2,1) * m(3,0) * m(3,2)
           -  m(3,1) * m(2,0) * m(3,2) + m(3,1) * m(3,0) * m(2,2);
  // cofactor (2,0) [2]
  inv(2,0) = m(1,0) * m(2,1) * m(3,3) - m(1,0) * m(3,1) * m(3,2)
           - m(1,1) * m(2,0) * m(3,3) + m(1,1) * m(3,0) * m(3,2)
           + m(3,1) * m(2,0) * m(3,1) - m(3,1) * m(3,0) * m(2,1);
  // cofactor (3,0) [3]
  inv(3,0) = -m(1,0) * m(2,1) * m(3,2) + m(1,0) * m(3,1) * m(2,2)
           +  m(1,1) * m(2,0) * m(3,2) - m(1,1) * m(3,0) * m(2,2)
           -  m(2,1) * m(2,0) * m(3,1) + m(2,1) * m(3,0) * m(2,1);
  // cofactor (1,1) [5]
  inv(1,1) = m(0,0) * m(2,2) * m(3,3) - m(0,0) * m(3,2) * m(3,2)
           - m(2,0) * m(2,0) * m(3,3) + m(2,0) * m(3,0) * m(3,2)
           + m(3,0) * m(2,0) * m(3,2) - m(3,0) * m(3,0) * m(2,2);
  // cofactor (2,1) [6]
  inv(2,1) = -m(0,0) * m(2,1) * m(3,3) + m(0,0) * m(3,1) * m(3,2)
           +  m(1,0) * m(2,0) * m(3,3) - m(1,0) * m(3,0) * m(3,2)
           -  m(3,0) * m(2,0) * m(3,1) + m(3,0) * m(3,0) * m(2,1);
  // cofactor (3,1) [7]
  inv(3,1) = m(0,0) * m(2,1) * m(3,2) - m(0,0) * m(3,1) * m(2,2)
           - m(1,0) * m(2,0) * m(3,2) + m(1,0) * m(3,0) * m(2,2)
           + m(2,0) * m(2,0) * m(3,1) - m(2,0) * m(3,0) * m(2,1);
  // cofactor (2,2) [10]
  inv(2,2) = m(0,0) * m(1,1) * m(3,3) - m(0,0) * m(3,1) * m(3,1)
           - m(1,0) * m(1,0) * m(3,3) + m(1,0) * m(3,0) * m(3,1)
           + m(3,0) * m(1,0) * m(3,1) - m(3,0) * m(3,0) * m(1,1);
  // cofactor (3,2) [11]
  inv(3,2) = -m(0,0) * m(1,1) * m(3,2) + m(0,0) * m(3,1) * m(2,1)
           +  m(1,0) * m(1,0) * m(3,2) - m(1,0) * m(3,0) * m(2,1)
           -  m(2,0) * m(1,0) * m(3,1) + m(2,0) * m(3,0) * m(1,1);
  // cofactor (3,3) [15]
  inv(3,3) = m(0,0) * m(1,1) * m(2,2) - m(0,0) * m(2,1) * m(2,1)
           - m(1,0) * m(1,0) * m(2,2) + m(1,0) * m(2,0) * m(2,1)
           + m(2,0) * m(1,0) * m(2,1) - m(2,0) * m(2,0) * m(1,1);

  // symetric
  inv(0,1) = inv(1,0); inv(0,2) = inv(2,0); inv(0,3) = inv(3,0);
  inv(1,2) = inv(2,1); inv(1,3) = inv(3,1);
  inv(2,3) = inv(3,2);

  // compute determinant
  Type det = m(0,0) * inv(0,0) + m(1,0) * inv(1,0) + m(2,0) * inv(2,0) + m(3,0) * inv(3,0);
  if (det == Type(0)) return Type(0);
  // inverse matrix
  inv /= det;
  return det;
}

/** @ingroup Algebra
 *  @brief compute the inverse of the symmetric 4x4 matrix m using only its lower part
 *  and store the result in inv.
 *  @note if the matrix is not invertible, the result will be a generalized inverse.
 *  @param m, inv the matrices to invert and its inverse
 *  @return The determinant value of m
 **/
template<class Matrix, int Size_>
static typename Matrix::Type invertSymMatrixXX( Matrix const& m, CArraySquare<typename Matrix::Type, Size_>& inv)
{
  typedef typename Matrix::Type Type;
#ifdef STKUSELAPACK
      lapack::SymEigen<Matrix> decomp(m);
      decomp.setUplo('L'); // default is U
#else
      SymEigen<Matrix> decomp(m);
#endif
  if (!decomp.run()) return Type(0);
  // compute tolerance
  Type tol = Arithmetic<Type>::epsilon() *decomp.norm();
  if (tol == 0) { tol = Arithmetic<Type>::min();}
  // compute (generalized) inverse matrix
  inv = decomp.rotation() * decomp.eigenValues().safeInverse(tol).diagonalize() * decomp.rotation().transpose();
  return decomp.det();
}


/** @ingroup hidden utility class allowing to call the correct static function
 *  computing the inverse of a symetric matrix.
 **/
template<class Matrix, int Size_>
struct invertSymMatrixDispatcher
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, Size_>& inv)
  { return invertSymMatrixXX(m, inv);}
};

template<class Matrix>
struct invertSymMatrixDispatcher<Matrix, 1>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 1>& inv)
  { return invertSymMatrix11(m, inv);}
};

template<class Matrix>
struct invertSymMatrixDispatcher<Matrix, 2>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 2>& inv)
  { return invertSymMatrix22(m, inv);}
};

template<class Matrix>
struct invertSymMatrixDispatcher<Matrix, 3>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 3>& inv)
  { return invertSymMatrix33(m, inv);}
};

template<class Matrix>
struct invertSymMatrixDispatcher<Matrix, 4>
{
  typedef typename Matrix::Type Type;
  inline static Type run( Matrix const& m, CArraySquare<Type, 4>& inv)
  { return invertSymMatrix44(m, inv);}
};

} // namespace hidden

} // namespace STK

#endif /* STK_SYMINVERT_H */
