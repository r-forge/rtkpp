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

/** @file STK_InvertFixedSizeMatrix.h
 *  @brief In this file we implement inversion method for fixed size matrix.
 **/

#ifndef STK_INVERTFIXEDSIZEMATRIX_H
#define STK_INVERTFIXEDSIZEMATRIX_H

#include <Arrays/include/STK_CArray.h>
#include "invert/STK_GenInvert.h"

namespace STK
{


/** @ingroup Algebra
 *  @brief The InvertSymMatrix class is a functor class allowing to compute the
 *  inverse of a symmetric matrix. It is specialized for fixed sized matrices.
 */
template<class Matrix, int Size> class InvertMatrix
{
  public:
    typedef typename Matrix::Type Type;
    /** Constructor */
    inline InvertMatrix( ArrayBase<Matrix> const& m)
                       : inv_(_R(0,Size-1))
                       , det_(hidden::invertMatrixDispatcher<Matrix, Size>::run(m, inv_))
    {}
    /** Destructor */
    inline virtual ~InvertMatrix() {}

    /** @return the inverse of the matrix m_ */
    inline CArraySquare<Type, Size> const& inv() const { return inv_;}
    /** @return the determinant of the matrix m_ */
    inline Type const& det() const { return det_;}
    /** @return @c true if the matrix m_ is invertible, @c false otherwise */
    inline bool isInvertible() const { return(det_!=0);}

    /** compute the inverse of the matrix m_. */
    inline CArraySquare<Type, Size> const& operator()() { return inv_;}

  protected:
    /** The inverse (or adjugate matrix if det_ is zero) of m_ */
    CArraySquare<Type, Size> inv_;
    /** determinant of the matrix m_ */
    Type det_;
};



}

#endif /* STK_INVERTFIXEDSIZEMATRIX_H_ */
