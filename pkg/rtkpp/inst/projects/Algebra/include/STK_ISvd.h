/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015  Serge Iovleff

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
 * Purpose:  Define The Interface ISvd Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_ISvd.h
 *  @brief In this file we define the interface class ISvd.
 **/
 
#ifndef STK_ISVD_H
#define STK_ISVD_H

#include "Sdk/include/STK_IRunner.h"
#include "Sdk/include/STK_IRecursiveTemplate.h"
#include "Arrays/include/STK_Array2D.h"
#include "Arrays/include/STK_Array2DSquare.h"
#include "Arrays/include/STK_Array2DDiagonal.h"

namespace STK
{
/** @ingroup Algebra
 *  @brief Compute the Singular Value Decomposition of an array.
 * 
 *  The method take as:
 *  - input: A matrix A(nrow,ncol)
 *  - output:
 *    -# U Array (nrow,ncolU).
 *    -# D Vector (ncol)
 *    -# V Array (ncol,ncol).
 *  and perform the decomposition: 
 *  - A = UDV'
 *  U can have more cols than A,
 *  and it is possible to compute some (all) vectors of Ker(A).
 **/
template<class Derived>
class ISvd  : public IRunnerBase, public IRecursiveTemplate<Derived>
{
  protected:
    typedef typename hidden::AlgebraTraits<Derived>::Array Array;
    /** Default constructor
     *  @param A the matrix to decompose.
     *  @param ref if true, U_ is a reference of A.
     *  @param withU if @c true save the left housolder transforms in @c U_.
     *  @param withV if @c true save the right housolder transforms in @c V_.
     **/
    ISvd( Array const& A, bool ref, bool withU = true, bool withV = true)
        : U_(A, ref), V_(), D_(), withU_(withU), withV_(withV), norm_(0), rank_(0)
    {}
    /** Copy Constructor
     *  @param S the Svd to copy
     **/
    ISvd( ISvd const& S)
        : U_(S.U_, S.U_.isRef()), V_(S.V_), D_(S.D_)
        , withU_(S.withU_), withV_(S.withV_)
        , norm_(S.norm_), rank_(S.rank_)
    {}
    /** destructor. */
    virtual ~ISvd() {}
    /** Operator = : overwrite the ISvd with S.
     *  @param S the Svd to copy
     **/
    ISvd& operator=(const ISvd &S)
    {
      U_ =S.U_; V_ = S.V_;  D_ = S.D_;
      withU_ =  S.withU_; withV_ = S.withV_;
      norm_ = S.norm_; rank_ = S.rank_;
      return *this;
    }
  public:
    /// @return the number of rows of U_
    inline int nrowU() const { return U_.sizeRows();}
    /// @return the number of columns of U_
    inline int ncolU() const { return U_.sizeCols();}
    /// @return the number of columns of D_
    inline int ncolD() const { return D_.sizeCols();}
    /// @return the number of rows of V_
    inline int nrowV() const { return V_.sizeRows();}
    /// @return the number of columns of V_
    inline int ncolV() const { return V_.sizeCols();}
    /// @return the norm of the matrix
    inline Real normSup()  const { return norm_;}
    /// @return the rank of the matrix
    inline int rank()  const { return rank_;}
    /// @return U
    inline Array const& getU() const { return U_;}
    /// @return  V
    inline ArraySquareX const& getV() const { return V_;}
    /// @return D
    inline ArrayDiagonalX const&  getD() const { return D_;}
    
  protected:
    /// U_ matrix
    Array U_;
    /// V_ square matrix
    ArraySquareX V_;
    /// Diagonal array of the singular values
    ArrayDiagonalX D_;
    /// Compute U_ ?
    bool withU_;
    /// Compute V_ ?
    bool withV_;
    /// norm of the matrix (largest singular value)
    Real norm_;
    /// rank of the matrix
    int  rank_;
};

} // namespace STK

#endif /* STK_ISVD_H */
