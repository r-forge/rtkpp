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
 * Project: stkpp::Arrays
 * created on: 07 jul. 2007
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SArray2D.h
  * @brief In this file, we define the final class @c SArray2D.
 **/
#ifndef STK_SARRAY2D_H
#define STK_SARRAY2D_H

#include "STK_SArray1D.h"

#include "STK_IArray2D.h"
#include "STK_IArray2DSlicers.h"
#include "STK_IArray2DModifiers.h"

namespace STK
{
// forward declaration
template< typename Type> class SArray2D;
template< typename Type> class SArray2DNumber;
template< typename Type> class SArray2DPoint;
template< typename Type> class SArray2DVector;

/** @ingroup Arrays
  * @brief Specialization of the SArray2D class for Real values.
  * An Array is a column oriented 2D container of Real.
 **/
typedef SArray2D<Real>   SArrayXX;
typedef SArray2D<double> SArrayXXd;
typedef SArray2D<int>    SArrayXXi;


namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for the SArray2D class.
 **/
template<class Type_>
struct Traits< SArray2D<Type_> >
{
  private:
    class Void {};
  public:
    typedef SArray2DNumber<Type_> Number;
    typedef SArray2DPoint<Type_>  Row;
    typedef SArray2DVector<Type_> Col;
    typedef SArray2DPoint<Type_>  SubRow;
    typedef SArray2DVector<Type_> SubCol;
    typedef SArray2D<Type_>       SubArray;
    typedef Void                 SubVector;

    typedef Type_                Type;
    typedef typename RemoveConst<Type>::Type const& TypeConst;

    enum
    {
      structure_ = Arrays::array2D_,
      orient_    = Arrays::by_col_,
      sizeRows_  = UnknownSize,
      sizeCols_  = UnknownSize,
      storage_   = Arrays::sparse_ // always dense
    };
    typedef SArray1D<Type, UnknownSize, UnknownSize> ColVector;
};

} // namespace hidden

/** @ingroup Arrays
  * @brief template sparse two dimensional column (vertically) oriented Array.
  *
  * A SArray2D is a sparse two-dimensional general implementation of an IArray2D.
  *
  * A column of an SArray2D is (almost like) an @c SArray2DVector
  * and a Row of an SArray2D is (almost like) an @c SArray2DPoint.
  *
  * When accessing to a row or a column, the methods return
  * the SArray2DPoint or SArray2DVector as a reference-like without copying
  * the data.
 *
 * The elements of the SArray2D can be acceded with the operators ()
 * - T(2,3) allow to access the third element of the second row of T
 * - T(3,Range(1,2)) allow to access to the first two members of the third
 *  row of T
 * - T(Range(2,3),2) allow to access to the second and third element of the
 *   second column of  T.
 * - T(Range(2,3),Range(1,2)) allow to access to the specified sub-matrix in T
 *
 * @sa SArray2DVector, SArray2DPoint, SArray2DSquare, SArray2DUpperTriangular, SArray2DLowerTriangular
 **/
template<class Type_ >
class SArray2D: public IArray2D< SArray2D<Type_> >
{
  public:
    /** Type for the Interface base Class. */
    typedef IArray2D< SArray2D<Type_> > Base;
    typedef ArrayBase < SArray2D<Type_> > LowBase;

    typedef typename hidden::Traits<SArray2D<Type_> >::Row Row;
    typedef typename hidden::Traits<SArray2D<Type_> >::Col Col;
    typedef typename hidden::Traits<SArray2D<Type_> >::SubRow SubRow;
    typedef typename hidden::Traits<SArray2D<Type_> >::SubCol SubCol;
    typedef typename hidden::Traits<SArray2D<Type_> >::SubVector SubVector;
    typedef typename hidden::Traits<SArray2D<Type_> >::SubArray SubArray;

    typedef typename hidden::Traits<SArray2D<Type_> >::Type Type;
    typedef typename hidden::Traits<SArray2D<Type_> >::TypeConst TypeConst;

    enum
    {
      structure_ = Arrays::array2D_,
      orient_    = Arrays::by_col_,
      sizeRows_  = UnknownSize,
      sizeCols_  = UnknownSize,
      storage_ = Arrays::dense_ // always dense
    };
    /** Default constructor */
    SArray2D(): Base() {}
    /** constructor
     *  @param I,J range of the rows and columns
     **/
    SArray2D( Range const& I, Range const& J): Base(I, J) {}
    /** constructor with rows and columns ranges specified and
     *  initialization with a constant.
     *  @param I,J range of the rows and columns
     *  @param v initial value of the container
     **/
    SArray2D( Range const& I, Range const& J, Type const& v): Base(I, J)
    { LowBase::setValue(v);}
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    SArray2D( SArray2D const& T, bool ref=false): Base(T, ref) {}
    /** Copy constructor by reference, ref_=1.
     *  @param T the container to wrap
     *  @param I,J range of the rows and columns to wrap
     **/
    template<class OtherArray>
    SArray2D( IArray2D<OtherArray> const& T, Range const& I, Range const& J)
          : Base(T, I, J) {}
    /** Copy constructor using an expression.
     *  @param T the container to wrap
     **/
    template<class OtherDerived>
    SArray2D( ExprBase<OtherDerived> const& T): Base()
    { LowBase::operator=(T);}
    /** Wrapper constructor Contruct a reference container.
     *  @param q pointer on the data
     *  @param I,J range of the rows and columns to wrap
     **/
    SArray2D( Type** q, Range const& I, Range const& J): Base(q, I, J) {}
    /** destructor. */
    ~SArray2D() {}
    /** operator = : overwrite the SArray2D with the right hand side T.
     *  @param T the container to copy
     **/
    template<class Rhs>
    SArray2D& operator=( ExprBase<Rhs> const& T) { return LowBase::operator=(T);}
    /** overwrite the SArray2D with T.
     *  @param T the container to copy
     **/
    SArray2D& operator=( SArray2D const& T) { return LowBase::assign(T);}
    /** set the container to a constant value.
     *  @param v the value to set
     **/
    SArray2D& operator=( Type const& v) { return LowBase::setValue(v);}
    /** New beginning index for the object.
     *  @param beg first index of the container
     **/
    void shift1D(int beg)
    { Base::shift(beg, beg);}
    /** New size for the container.
     *  @param I range of the columns and rows of the container
     **/
    SArray2D& resize1D( Range const& I)
    { Base::resize(I, I); return *this;}
    /** Swapping the pos1 row and the pos2 row.
     *  @param pos1,pos2 position of the rows to swap
     **/
    void swapRows( int const& pos1, int const& pos2)
    {
#ifdef STK_BOUNDS_CHECK
      // check conditions
      if (this->beginRows() > pos1)
      STKOUT_OF_RANGE_2ARG(SArray2D::swapRows,pos1, pos2,beginRows()>pos1);
      if (this->endRows() <= pos1)
      STKOUT_OF_RANGE_2ARG(SArray2D::swapRows,pos1, pos2,endRows()<=pos1);
      if (this->beginRows() > pos2)
      STKOUT_OF_RANGE_2ARG(SArray2D::swapRows,pos1, pos2,beginRows()>pos2);
      if (this->endRows() <= pos2)
      STKOUT_OF_RANGE_2ARG(SArray2D::swapRows,pos1, pos2,endRows()<=pos2);
#endif
      for (int j=this->beginCols(); j<this->endCols(); j++)
      { std::swap(this->elt(pos1, j), this->elt(pos2, j));}
    }
};

} // namespace STK

#endif
// STK_SARRAY2D_H