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

/** @file STK_InvertSymMatrix.h
 *  @brief In this file we implement inversion method for fixed size symmetric matrices.
 **/

#ifndef STK_INVERTSYMMATRIX_H
#define STK_INVERTSYMMATRIX_H

#include <Arrays/include/STK_CArraySquare.h>
#include "invert/STK_SymInvert.h"

namespace STK
{

/** @ingroup Algebra
 *  @brief The InvertSymMatrix class is a functor class allowing to compute the
 *  inverse of a symmetric matrix. It is specialized for fixed sized matrices
 *  of small size (up to size=4).
 *  For general matrix, it uses the eigenvalues decomposition of symmetric matrix
 *  and compute the inverse using the formula \f$ PD^{-1}P^{-1} \f$.
 */
template<class Matrix, int Size_>
class InvertSymMatrix
{
  public:
    typedef typename Matrix::Type Type;
    /** Constructor */
    inline InvertSymMatrix( ArrayBase<Matrix> const& m)
                          : inv_(_R(0,Size_-1))
                          , det_(hidden::invertSymMatrixDispatcher<Matrix, Size_>::run(m.asDerived(), inv_))
    { inv_.shift(m.beginRows());}
    /** Destructor */
    inline ~InvertSymMatrix() {}

    /** @return the inverse of the matrix m_ */
    inline CArraySquare<Type, Size_> const& inv() const { return inv_;}
    /** @return the determinant of the matrix m_ */
    inline Type const& det() const { return det_;}
    /** @return @c true if the matrix m_ is invertible, @c false otherwise */
    inline bool isInvertible() const { return det_!=0;}

    /** Operator().
     *  @return the inverse of the matrix m_. */
    inline CArraySquare<Type, Size_> const& operator()() { return inv_;}

  protected:
    /** The inverse (or adjugate matrix if det_ is zero) of m_ */
    CArraySquare<Type, Size_> inv_;
    /** determinant of the matrix m_ */
    Type det_;
};

/** @ingroup Algebra
 *  @brief The GInvertSymMatrix class is a functor allowing to compute the generalized
 *  inverse of a symmetric matrix.
 */
class GInvertSymMatrix
{
  public:
    /** Constructor */
    inline GInvertSymMatrix() {}
    /** Destructor */
    inline ~GInvertSymMatrix() {}
    /** compute the generalized inverse of the symmetric matrix x. x is
     *  overwritten by the result.
     *  @param x the matrix to inverse.
     **/
    template<class Matrix>
    inline Matrix& operator()(ArrayBase<Matrix>*& x)
    {
      enum
      {
        sizeRows_ = Matrix::sizeRows_,
        sizeCols_ = Matrix::sizeCols_,
        size_ = (sizeRows_ < sizeCols_) ? sizeRows_ : sizeCols_
      };
      InvertSymMatrix<Matrix, size_> funct(x->asDerived());
      if ((funct.det()==0)&&(size_ != UnknownSize))
      {
#ifdef STKUSELAPACK
        lapack::SymEigen<Matrix> decomp(x->asDerived());
#else
        SymEigen<Matrix> decomp(x->asDerived());
#endif
        decomp.run();
        decomp.ginv(x->asDerived());
      }
      else
      { x->asDerived() = funct.inv();}
      return x->asDerived();
    }
    /** Compute the generalized inverse of the symmetric matrix x. x is
     *  overwritten by the result.
     *  @param x the matrix to inverse.
     **/
    template<class Matrix>
    inline Matrix& operator()(ArrayBase<Matrix>& x)
    {
      enum
      {
        sizeRows_ = Matrix::sizeRows_,
        sizeCols_ = Matrix::sizeCols_,
        size_ = (sizeRows_ < sizeCols_) ? sizeRows_ : sizeCols_
      };
      InvertSymMatrix<Matrix, size_> funct(x.asDerived());
      if ((funct.det()==0)&&(size_ != UnknownSize))
      {
#ifdef STKUSELAPACK
        lapack::SymEigen<Matrix> decomp(x.asDerived());
#else
        SymEigen<Matrix> decomp(x.asDerived());
#endif
        decomp.run();
        decomp.ginv(x.asDerived());
      }
      else
      { x.asDerived() = funct.inv();}
      return x.asDerived();
    }
};

} // namespace STK

#endif /* STK_INVERTFIXEDSIZESYMMATRIX_H */
