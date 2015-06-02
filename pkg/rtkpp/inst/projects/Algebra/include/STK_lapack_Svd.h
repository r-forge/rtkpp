/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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
 * created on: 20 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_lapack_Svd.h
 *  @brief In this file we define the enclosing class of the dgeqrf lapack routine.
 **/


#ifndef STK_LAPACK_SVD_H
#define STK_LAPACK_SVD_H

//#include "STK_ISvd.h"

#ifdef STKUSELAPACK

extern "C"
{
#ifdef STKREALAREFLOAT
/** LAPACK routine in float to compute the SVD decomposition */
extern void sgeqrf_(int* M, int* N, float* A, int* LDA, float* TAU, float* WORK, int* LWORK, int* INFO );
#else
/** LAPACK routine in double to compute the SVD decomposition */
extern void dgeqrf_(int* M, int* N, double* A, int* LDA, double* TAU, double* WORK, int* LWORK, int* INFO );
#endif
}

#endif // STKUSELAPACK

namespace STK
{

namespace lapack
{
/** @ingroup Algebra
 *  {
 *    @class Svd
 *    @brief Svd computes the SVD decomposition of a real matrix using the
 *    Lapack routine dgeqrf.
 */
class Svd : public ISvd<Svd>
{
  public:
    typedef ISvd<Svd> Base;
    /** Default constructor.
     *  @param data the matrix to decompose
     *  @param ref true if we overwrite A
     **/
    inline Svd( ArrayXX const&  data, bool ref = false): Base(data, ref) {}
    /** @brief Constructor
     *  @param data reference on a matrix expression
     */
    template<class Derived>
    Svd( ArrayBase<Derived> const& data): Base(data){}
    /** Copy constructor.
     *  @param decomp the decomposition  to copy
     **/
    inline Svd( Svd const& decomp): Base(decomp) {}
    /** virtual destructor */
    inline virtual ~Svd() {}
    /** @brief clone pattern */
    inline virtual Svd* clone() const { return new Svd(*this);}
    /** Operator = : overwrite the Svd with decomp. */
    inline Svd& operator=(Svd const& decomp)
    {
      Base::operator=(decomp);
      return *this;
    }
    /** @brief Run qr decomposition
     *  Launch geqrf LAPACK routine to perform the qr decomposition.
     *  @return @c true if no error occur, @c false otherwise
     */
    bool runImpl();

  protected:
    /** wrapper of the LAPACK DGESVDF routine. Compute the Svd decomposition
     *  of a matrix.
     *
     * @param[in] m The number of rows of the matrix A.  M >= 0.
     *
     * @param[in] n The number of columns of the matrix A.  N >= 0.
     *
     * @param[in,out] a Real array, dimension (LDA, N)
     * \verbatim
     *     On entry, the M-by-N matrix A.
     *     On exit, the elements on and above the diagonal of the array
     *     contain the min(M,N)-by-N upper trapezoidal matrix R (R is
     *     upper triangular if m >= n); the elements below the diagonal,
     *     with the array TAU, represent the orthogonal matrix Q as a
     *     product of min(m,n) elementary reflectors (see Further Details).
     * \endverbatim
     *
     * @param[in] lda The leading dimension of the array A.  LDA >= max(1,M).
     *
     * @param[out] tau Real array, dimension min(M,N)
     * The scalar factors of the elementary reflectors (see Further Details).
     *
     * @param[in,out] work Real array, dimension (MAX(1,LWORK))
     * \verbatim
     *   On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
     * \endverbatim
     *
     * @param[in] lwork The  dimension  of  the array WORK
     * \verbatim
     *  LWORK >= max(1,N).
     *  For optimum performance LWORK >= N*NB, where NB is the optimal blocksize.
     *
     *  If LWORK = -1, then a workspace query is assumed; the routine
     *  only calculates the optimal size of the WORK array, returns
     *  this value as the first entry of the WORK array, and no error
     *  message related to LWORK is issued by XERBLA.
     * \endverbatim
     *
     * @return info
     * \verbatim
     *  = 0:  successful exit
     *  < 0:  if INFO = -i, the i-th argument had an illegal value
     * \endverbatim
     *
     * @verbatim
     *  Further Details
     *  ===============
     *
     *  The matrix Q is represented as a product of elementary reflectors
     *
     *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
     *
     *  Each H(i) has the form
     *
     *     H(i) = I - tau * v * v'
     *
     *  where tau is a real scalar, and v is a real vector with
     *  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
     *  and tau in TAU(i).
     * @endverbatim
     */
    int geqrf(int m, int n, Real* a, int lda, Real* tau, Real *work, int lwork);
};


/** @} */

} // namespace lapack

} // namespace STK

#endif /* STK_LAPACK_SVD_H */
