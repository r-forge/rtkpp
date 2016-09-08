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

/** @file STK_InvertFixedSizeSymMatrix.h
 *  @brief In this file we implement inversion method for fixed size symmetric matrices.
 **/

#ifndef STK_INVERTFIXEDSIZESYMMATRIX_H
#define STK_INVERTFIXEDSIZESYMMATRIX_H

#include <Arrays/include/STK_CArray.h>

namespace STK
{

/** @ingroup Algebra
 *  @brief compute the inverse of the symmetric 1x1 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix11( Matrix const& m, CArraySquare<typename Matrix::Type, 1>& inv);
/** @ingroup Algebra
 *  @brief compute the inverse of the symmetric 2x2 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix22( Matrix const& m, CArraySquare<typename Matrix::Type, 2>& inv);
/** @ingroup Algebra
 *  @brief compute the inverse of the symmetric 3x3 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix33( Matrix const& m, CArraySquare<typename Matrix::Type, 3>& inv);
/** @ingroup Algebra
 *  @brief compute the inverse of the symmetric 4x4 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
 **/
template<class Matrix>
static typename Matrix::Type invertSymMatrix44( Matrix const& m, CArraySquare<typename Matrix::Type, 4>& inv);

namespace hidden
{
/** @ingroup hidden utility class allowing to call the correct static function
 *  computing the inverse of a symetric matrix. */
template<class Matrix, int Size> struct invertSymMatrixDispatcher;

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

/* compute the inverse of the symmetric 2x2 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
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
/* compute the inverse of the symmetric 2x2 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
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
/** compute the inverse of the symmetric 3x3 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
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
/* compute the inverse of the symmetric 4x4 matrix m using only its lower part
 *  and store the result in inv.
 *  @param m, inv the matrices to invert and its inverse
 *  @return @c true if m is invertible, @c false otherwise
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
 *  @brief The InvertSymMatrix class is a functor class allowing to compute the
 *  inverse of a symmetric matrix. It is specialized for fixed sized matrices.
 */
template<class Matrix, int Size> class InvertSymMatrix
{
  public:
    typedef typename Matrix::Type Type;
    /** Constructor */
    inline InvertSymMatrix( ArrayBase<Matrix> const& m)
                          : m_(m.asDerived())
                          , inv_(_R(0,Size-1))
                          , det_(hidden::invertSymMatrixDispatcher<Matrix, Size>::run(m_, inv_))
                          , isInvertible_(det_!=0)
    {}
    /** Destructor */
    inline virtual ~InvertSymMatrix() {}

    /** @return the inverse of the matrix m_ */
    inline CArraySquare<Type, Size> const& inv() const { return inv_;}
    /** @return the determinant of the matrix m_ */
    inline Type const& det() const { return det_;}
    /** @return @c true if the matrix m_ is invertible, @c false otherwise */
    inline bool const& isInvertible() const { return isInvertible_;}

    /** compute the inverse of the matrix m_. */
    inline CArraySquare<Type, Size> const& operator()() { return inv_;}

  protected:
    /** a constant reference on the 4x4 matrix m_ to invert */
    Matrix const& m_;
    /** The inverse (or adjugate matrix if det_ is zero) of m_ */
    CArraySquare<Type, Size> inv_;
    /** determinant of the matrix m_ */
    Type det_;
    /** @c true if the matrix m_ is invertible */
    bool isInvertible_;
};

} // namespace STK

#endif /* STK_INVERTFIXEDSIZESYMMATRIX_H */
