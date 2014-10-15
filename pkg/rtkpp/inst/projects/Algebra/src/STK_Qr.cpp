/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2007  Serge Iovleff

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
 * Project:  Algebra
 * Purpose:  Implement the Qr Class.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 *
 **/

/** @file STK_Qr.cpp
 *  @brief In this file we implement the Qr Class (QR decomposition).
 **/
 
#include "../include/STK_Qr.h"

#include "Arrays/include/STK_Array2DPoint.h"
//#include "Arrays/include/STK_Array2D_Functors.h"

#include "../include/STK_Householder.h"
#include "../include/STK_Givens.h"

#ifdef STK_ALGEBRA_DEBUG
#include "../../Arrays/include/STK_Display.h"

template< class Container2D >
void print(Container2D const& A, STK::String const& name)
{
  stk_cout << "print: " << name << _T("\n";);
  stk_cout << name << _T(".isRef() =")        << A.isRef()  << _T("\n");
  stk_cout << name << _T(".capacityHo() =")   << A.capacityHo()  << _T("\n");
  stk_cout << name << _T(".cols() =")      << A.cols()  << _T("\n");
  stk_cout << name << _T(".rows() =")      << A.rows()  << _T("\n\n");
  stk_cout << name << _T(".rangeCols().isRef() =")  << A.rangeCols().isRef() << _T("\n");
  stk_cout << name << _T(".rangeCols() =\n")  << A.rangeCols() << _T("\n");
  stk_cout << name << _T(".capacityCols().isRef() =") << A.capacityCols().isRef()  << _T("\n");
  stk_cout << name << _T(".capacityCols() =\n") << A.capacityCols()  << _T("\n");
  stk_cout << name << _T(".p_data() =") << A.allocator().p_data()  << _T("\n");
}
#endif

namespace STK
{

/* Constructor */
Qr::Qr( Matrix const& A, bool ref) : Q_(A, ref), R_(), compq_(false)
{ run();}

/* Computing the QR decomposition of the matrix Q_. */
void Qr::run()
{
  if (Q_.empty())     // if the container is empty
  {
    compq_ = true;  // Q_ is computed
    return;
  }
  // translate the beg to 1
  Q_.shift(1,1);
  qr();
}


/* Computation of the QR decomposition */
void Qr::qr()
{
  R_.resize(Q_.rows(), Q_.cols());
  R_ = 0.0;
  // start left householder reflections
  Range r(Q_.rows()), c(Q_.cols());
  for(int j = R_.beginRows(); (j < R_.endRows()) && (j < Q_.endCols()) ; ++j)
  {
    Vector u(Q_, r, j);    // get a reference on the j-th column in the range r
    R_(j, j) = house(u);   // compute the housolder vector
    leftHouseholder(Q_(r, c.incFirst(1)), u); // apply-it to the remaining part of Q_
    r.incFirst(1);
    R_.row(j, c).copy(Q_.row(j, c)); // copy current row in R_
  }
}


/* Computation of Q. */
void Qr::compQ()
{
#ifdef STK_ALGEBRA_VERBOSE
  stk_cout << _T("Entering Qr::compQ(). compq_ =") << compq_ << _T("\n");
#endif
  // if Q_ is computed yet
  if (compq_) return;
  // number of non zero cols of Q_  
  int ncol  = std::min(Q_.sizeRows(), Q_.sizeCols()), lastCol;
  // add or remove the column
  if (ncol < Q_.sizeCols())
  {
    Q_.popBackCols(Q_.sizeCols() - ncol);
    lastCol = Q_.lastIdxCols();
  }
  else
  {
    lastCol = Q_.lastIdxCols();
    stk_cout << "ncol = " << ncol << " lastCol = " << lastCol << "\n";
    if (ncol < Q_.sizeRows())
    {
      Q_.pushBackCols(Q_.sizeRows() -ncol);
      // Initialize added columns
      Q_.col( _R( lastCol+1, Q_.lastIdxCols()) ).setValue(0);
      for (int i=lastCol+1; i< Q_.endCols(); ++i) { Q_(i, i) = 1.0;}
    }
  }
  // compute other columns
  for (int iter=lastCol, iter1= lastCol +1; iter>=Q_.beginCols(); iter--, iter1--)
  {
    // get current householder vector
    Vector u(Q_, _R(iter, Q_.lastIdxRows()), iter);
    // Apply Householder vector to the right of the matrix
    leftHouseholder( Q_( _R(iter, Q_.lastIdxRows()), _R(iter1, Q_.lastIdxCols())), u);
    // Get beta and test
    Real beta = Q_(iter,iter);
    // update the column iter
    Q_(iter,iter) = 1.0 + beta;
    Q_(_R(iter1, Q_.lastIdxRows()), iter ) *= beta;
    // update the column iter
    Q_( Range(1,iter-1, 0), iter) = 0.0;
  }
  // Q_ is now computed
  compq_ = true;
#ifdef  STK_ALGEBRA_VERBOSE
  stk_cout << _T("Terminating Qr::compQ(). compq_ =") << compq_ << _T("\n");
#endif
}


/* Destructor */
Qr::~Qr() {}

/* clear Q_ and R_. */
void Qr::clear()
{
  Q_.clear();
  R_.clear();
}


/* Operator = : overwrite the Qr with S. */
Qr& Qr::operator=(const Qr& S)
{
  compq_ = S.compq_;       //< Is Q computed ?
  Q_ = S.Q_;               //< Matrix V
  R_ = S.R_;               //< Singular values

  return *this;
}

/* Delete the jth column and update the QR decomposition : default
 * is the last col
 **/    
void Qr::popBackCols(int const& n)
{
  // delete n cols
  R_.popBackCols(n);
}

void Qr::eraseCol(int const& pos)
{
  if (pos < R_.beginCols())
  { STKOUT_OF_RANGE_1ARG(Qr::eraseCol,pos,pos<R_.beginCols());}
  if (R_.lastIdxCols() < pos)
  { STKOUT_OF_RANGE_1ARG(Qr::eraseCol,pos,pos<R_.lastIdxCols()<pos);}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // compute the number of iteration for updating to zeroed
  int niter = R_.beginCols()-1+std::min(R_.sizeRows(), R_.sizeCols());
  // Zeroed the unwanted elements (z)
  for (int iter = pos+1; iter<=niter; iter++)
  {
    Real sinus, cosinus;
    // compute the Givens rotation
    R_(iter-1, iter) = compGivens( R_(iter-1, iter), R_(iter, iter), cosinus, sinus);
    R_(iter, iter)   = 0.0;
    // if necessary update R_ and Q_
    if (sinus)
    {
      // create a reference on the sub-Matrix
      MatrixUpperTriangular Rsub(R_.col(Range(iter+1, R_.lastIdxCols(), 0)), true);
      // Update the next rows (iter1:ncolr_) of R_
      leftGivens (Rsub, iter-1, iter, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iter-1, iter, cosinus, sinus);
    }
  }
  // erase the column pos
  R_.eraseCols(pos);
  
  // update the range of the remaining cols of the container
  R_.update(Range(pos, std::min(R_.lastIdxRows(), R_.lastIdxCols()), 0));
}


/* Adding the last column and update the QR decomposition. */
void Qr::pushBackCol(Vector const& T)
{
  // check conditions
  if (T.range() != Q_.rows())
  { STKRUNTIME_ERROR_NO_ARG(Qr::pushBackCol,T.range() != Q_.rows());}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  int lastColR = R_.endCols();
  // Create an auxiliary container
  Vector Rncolr = Q_.transpose() * T; // Rncolr of size Q_.cols()
  // update Q_
  for (int iter = Q_.lastIdxCols()-1, iter1 = Q_.lastIdxCols(); iter>=lastColR; iter--, iter1--)
  { 
    Real sinus, cosinus;
    // compute the Givens rotation
    Rncolr[iter] = compGivens( Rncolr[iter], Rncolr[iter1], cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    { rightGivens(Q_, iter, iter1, cosinus, sinus);}
  }
  // update R_
  R_.pushBackCols();
  R_.col(lastColR).copy(Rncolr.sub(R_.rangeRowsInCol(lastColR)));
}


/* Adding the jth column and update the QR decomposition.   */
void Qr::insertCol(Vector const& T, int const& pos)
{
  if (pos < R_.beginCols())
  { STKOUT_OF_RANGE_1ARG(Qr::insertCol,pos,pos<R_.beginCols());}
  if (R_.lastIdxCols() < pos)
  { STKOUT_OF_RANGE_1ARG(Qr::insertCol,pos,pos<R_.lastIdxCols()<pos);}
  if (T.range() != Q_.rows())
  { STKRUNTIME_ERROR_1ARG(Qr::insertCol,pos,T.range() != Q_.rows());}
  // if Q_ is not computed yet
  if (!compq_) compQ();
  // Adding a column to R
  R_.insertCols(pos);
  // update the range of the remaining cols of R_
  R_.update( _R(pos+1, std::min(R_.lastIdxRows(), R_.lastIdxCols())) );
  for (int i=pos+1; i< std::min(R_.endRows(), R_.endCols()); ++i) R_(i,i) = 0.0;

  Vector Rpos =  Q_.transpose() * T;
  // Zeroed the unwanted elements
  for (int iter= Q_.lastIdxCols(), iterm1= Q_.lastIdxCols()-1; iter>pos; iterm1--, iter--)
  { 
    Real sinus, cosinus;
    // compute the Givens rotation
    Rpos[iterm1]  = compGivens(Rpos[iterm1], Rpos[iter], cosinus, sinus);
    // apply Givens rotation if necessary
    if (sinus)
    {
      // create a reference on the sub-Matrix
      MatrixUpperTriangular Rsub(R_.col(_R(iter, R_.lastIdxCols())), true);
      // Update the next rows (iter:ncolr_) of R_
      leftGivens( Rsub, iterm1, iter, cosinus, sinus);
      // Update the cols of Q_
      rightGivens(Q_, iterm1, iter, cosinus, sinus);
    }
  }
  // update R_
  R_.col(pos) = Rpos.sub(R_.rangeRowsInCol(pos));
  R_.update(pos);
}

} // namespace STK

