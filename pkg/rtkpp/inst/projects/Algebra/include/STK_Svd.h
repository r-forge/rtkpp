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
 * Purpose:  Define The Svd Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Svd.h
 *  @brief In this file we define the Svd Class.
 **/
 
#ifndef STK_SVD_H
#define STK_SVD_H

#include "STK_ISvd.h"
#include "Arrays/include/STK_Array2DPoint.h"
#include "Arrays/include/STK_Array2DVector.h"

namespace STK
{
// forward declaration
//template<class Array> class Svd;
//
//namespace hidden
//{
///** @ingroup hidden
// *  Specialization for the Svd class.
// **/
//template<class Array_>
//struct AlgebraTraits< Svd<Array_> >
//{
//  typedef Array_ Array;
//};
//
//} // namespace hidden

/** @ingroup Algebra
 *  @brief The class Svd compute the Singular Value Decomposition
 *  of a Array with the Golub-Reinsch Algorithm.
 * 
 *  The method take as:
 *  - input: a matrix A(nrow,ncol)
 *  - output:
 *    -# U Array (nrow,ncol).
 *    -# D diagonal matrix (min(norw,ncol))
 *    -# V Array (ncol,ncol).
 *  and perform the decomposition: 
 *  - A = UDV' (transpose V).
 *  U can have more columns than A,
 *  and it is possible to ompute some (all) vectors of Ker(A).
 **/
class Svd : public ISvd<Svd>
{
  public :   
    typedef ISvd<Svd> Base;
    /** Default constructor
     *  @param A the matrix to decompose.
     *  @param ref if true, U_ is a reference of A.
     *  @param withU if true, we save the left housolder transforms in U_.
     *  @param withV if true, we save the right housolder transforms in V_.
     **/
    Svd( ArrayXX const& A, bool ref= false, bool withU= true, bool withV= true);
    /** Copy Constructor
     *  @param S the Svd to copy
     **/
    Svd( const Svd &S);
    /** destructor. */
    virtual ~Svd();
    /** Operator = : overwrite the Svd with S.
     *  @param S the Svd to copy
     **/
    Svd& operator=(const Svd &S);
    /** run the Svd */
    virtual bool run();
    /** Compute the svd of the Array A and copy the data
     *  see the corresponding constructor Take care that if U_ was previously
     *  a reference, it cannot be modified.
     *  @param A is the matrix to decompose.
     *  @param withU if true, we save the left housolder transforms
     *  in U_.
     *  @param withV if true, we save the right housolder transforms
     *  in V_.
     **/    
    void setData( ArrayXX const& A, bool withU = true, bool withV = true);
    /** Computing the bidiagonalization of M.
     *  The diagonal and the subdiagonal are stored in D and F
     *  @param M the matrix to bi-diagonalize, the matrix is overwritten
     *  with the left and right Householder vectors.
     *  The method return the estimate of the inf norm of M.
     *  @param D the element of the diagonal
     *  @param F the element of the surdiagnal
     **/
    static Real bidiag(const ArrayXX& M, ArrayDiagonalX& D, Vector& F);
    /** right eliminate the element on the subdiagonal of the row nrow
     *  @param D the diagonal of the matrix
     *  @param F the subdiagonal of the matrix
     *  @param nrow the number of the row were we want to rightEliminate
     *  @param V a right orthogonal Square ArrayXX
     *  @param withV true if we want to update V
     *  @param tol the tolerance to use
     **/
    static void rightEliminate( ArrayDiagonalX& D
                              , Vector& F
                              , int const& nrow
                              , ArraySquareX& V
                              , bool withV = true
                              , Real const& tol = Arithmetic<Real>::epsilon()
                              );
    /** left eliminate the element on the subdiagonal of the row @c nrow
     *  @param D the diagonal of the matrix
     *  @param F the subdiagonal of the matrix
     *  @param nrow the number of the row were we want to rightEliminate
     *  @param U a left orthogonal ArrayXX
     *  @param withU true if we want to update U
     *  @param tol the tolerance to use
     **/
    static void leftEliminate( ArrayDiagonalX& D
                             , Vector& F
                             , int const& nrow
                             , ArrayXX& U
                             , bool withU = true
                             , Real const& tol = Arithmetic<Real>::epsilon()
                             );

    /** Computing the diagonalization of a bi-diagonal matrix
     *  @param D the diagonal of the matrix
     *  @param F the subdiagonal of the matrix
     *  @param U a left orthogonal ArrayXX
     *  @param withU true if we want to update U
     *  @param V a right orthogonal Square ArrayXX
     *  @param withV true if we want to update V
     *  @param tol the tolerance to use
     **/
    static bool diag( ArrayDiagonalX& D
                    , Vector& F
                    , ArrayXX& U
                    , ArraySquareX& V
                    , bool withU = true
                    , bool withV = true
                    , Real const& tol = Arithmetic<Real>::epsilon()
                    );
  private:
    /// Values of the Sub-diagonal
    Vector F_;
    /// Initialize the containers
    void init();
    /// Svd main steps
    bool compSvd();
    /// Compute U (if withU_ is true)
    void compU();
    /// Compute V (if withV_ is true)
    void compV();
};

} // namespace STK

#endif
// STK_SVD_H
