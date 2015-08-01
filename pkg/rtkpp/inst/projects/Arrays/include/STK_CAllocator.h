/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015 Serge Iovleff

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
 * Project:  stkpp::Arrays
 * created on: 10 août 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_CAllocator.h
 *  @brief In this file we define the CAllocator class.
 **/

#ifndef STK_CALLOCATOR_H
#define STK_CALLOCATOR_H

#include "STK_ITContainer2D.h"
#include "STK_AllocatorBase.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief Allocator for dense Array classes.
 *  The data are stored in two dimensions.
 *  It can be the columns or the rows allocator of any dense container.
 */
template<typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CAllocator.
 */
template< typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
struct Traits< CAllocator<Type_, SizeRows_, SizeCols_, Orient_> >
{
  private:
    class Void { };

  public:
    enum
    {
      orient_    = Orient_,
      sizeRows_  = SizeRows_,
      sizeCols_  = SizeCols_,
      sizeProd_  = (SizeRows_ >= SqrtUnknownSize)||(SizeCols_ >= SqrtUnknownSize) ? UnknownSize : SizeRows_ * SizeCols_,
      storage_   = Arrays::dense_ // always dense
    };

    typedef CAllocator<Type_, 1, 1, Orient_> Number;

    typedef CAllocator<Type_, 1, SizeCols_, Orient_> Row;
    typedef CAllocator<Type_, SizeRows_, 1, Orient_> Col;

    typedef CAllocator<Type_, 1, UnknownSize, Orient_> SubRow;
    typedef CAllocator<Type_, UnknownSize, 1, Orient_> SubCol;

    /** If one of the Size is 1, we have a Vector (a column) or a Point (a row)
     *  (What to do if both are = 1 : Type_ or array (1,1) ?).
     **/
    typedef typename If< (SizeRows_ == 1)||(SizeCols_ == 1)   // one row or one column
                       , typename If<(SizeCols_ == 1)
                                    , typename If<SizeRows_==1, Number, SubCol>::Result
                                    , SubRow>::Result
                       , Void
                       >::Result SubVector;
    /** If one of the Size is 1, we have a Vector (a column) or a Point (a row)
     *  (What to do if both are = 1 : Type_ or array (1,1) ?).
     **/
    typedef typename If< (SizeRows_ >= SqrtUnknownSize)||(SizeCols_ >= SqrtUnknownSize)   // one row or one column
                       , AllocatorBase<Type_, UnknownSize>
                       , AllocatorBase<Type_, sizeProd_>
                       >::Result Allocator;
    typedef Type_  Type;
    typedef typename RemoveConst<Type_>::Type const& ReturnType;

    // use this as default. FIXME: Not optimal in case we just get a SubArray
    // with unmodified rows or cols size.
    typedef CAllocator<Type, UnknownSize, UnknownSize, Orient_> SubArray;
};

} // namespace hidden

// forward declaration
template<class Derived, bool Orient_> class OrientedCAllocator;
template<class Derived, int SizeRows_, int SizeCols_, bool Orient_> class StructuredCAllocator;

/**  @ingroup Arrays
 *   @brief Specialization for column-oriented Allocators.*/
template<class Derived>
class OrientedCAllocator<Derived, Arrays::by_col_>: public ITContainer2D<Derived>
{
  protected:
    typedef ITContainer2D<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;

    /** default constructor */
    inline OrientedCAllocator( Range const& I, Range const& J)
                             : Base(I, J), ldx_(I.size()), allocator_(prod(I, J))
    {}
    /** copy constructor */
    inline OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                             : Base(A), ldx_(A.ldx()), allocator_(A.allocator_, ref)
    { if (!ref) allocator_.copy(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, Arrays::by_col_> const& A
                             , Range const& I, Range const& J)
                             : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)  {}
    /** wrapper constructor for 0 based C-Array*/
    inline OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                             : Base(Range(0,nbRow), Range(0,nbCol)), ldx_(nbRow)
                             , allocator_(q, nbRow * nbCol, true)
    {}
    /** destructor */
    inline ~OrientedCAllocator() {}

  public:
    Type* p_data() { return allocator_.p_data();}
    Type const* p_data() const { return allocator_.p_data();}
    bool isRef() const { return allocator_.isRef();}
    Allocator& allocator() { return allocator_;}
    Allocator const& allocator() const { return allocator_;}

    /** @return the index of the allocator*/
    inline int ldx() const { return ldx_;}
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline Type const& elt2Impl(int i, int j) const { return p_data()[j*ldx_ + i];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i, j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return p_data()[j*ldx_ + i];}
    /** set a value to this allocator.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int j= this->beginCols(); j < this->endCols(); ++j)
        for (int i = this->beginRows(); i < this->endRows(); ++i)
        { this->elt(i, j) = v;}
    }
  protected:
    /** index of the data set */
    int ldx_;
    /** set index of the data. */
    inline void setIdx( int idx) { ldx_ = idx;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(OrientedCAllocator &T)
    {
      allocator_.exchange(T.allocator_);
      Base::exchange(T);
      std::swap(ldx_, T.ldx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with  range I for the rows and range J for the
     *  columns.
     *  @param I,J the range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static inline Range prod(Range const& I, Range const& J)
    { return Range(I.size()*J.begin()+I.begin(), I.size()*J.size()); }
    /** return the increment to apply to a zero based pointer corresponding to
     *  the actual first row and first column indexes. */
    inline int shiftInc(int firstRow, int firstCol)
    { return ldx_*firstCol+firstRow; }
    /** set the index corresponding to the actual size of the allocator. */
    inline void setSizedIdx() {ldx_ = this->asDerived().sizeRows();}
    /** manager of the memory */
    Allocator allocator_;
};

/**  @ingroup Arrays
 *   @brief Specialization for row-oriented Allocators.*/
template<class Derived>
class OrientedCAllocator<Derived, Arrays::by_row_>: public ITContainer2D<Derived>
{
  protected:
    typedef ITContainer2D<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Allocator Allocator;
    /** constructor with specified ranges */
    inline OrientedCAllocator( Range const& I, Range const& J)
                             : Base(I, J), ldx_(J.size()), allocator_(prod(I, J))
    {}
    /** copy constructor */
    inline OrientedCAllocator( OrientedCAllocator const& A, bool ref)
                             : Base(A), ldx_(A.ldx_), allocator_(A.allocator_, ref)
    { if (!ref) allocator_.copy(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedCAllocator( OrientedCAllocator<OtherDerived, Arrays::by_row_> const& A
                             , Range const& I, Range const& J)
                             : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline OrientedCAllocator( Type* const& q, int nbRow, int nbCol)
                             : Base(Range(0,nbRow), Range(0,nbCol)), ldx_(nbCol)
                             , allocator_(q, nbRow * nbCol, true)
    {}
    /** destructor */
    inline ~OrientedCAllocator() {}

  public:
    Type* p_data() { return allocator_.p_data();}
    Type const* p_data() const { return allocator_.p_data();}
    bool isRef() const { return allocator_.isRef();}
    Allocator& allocator() { return allocator_;}
    Allocator const& allocator() const { return allocator_;}

    /** @return the index of the allocator*/
    inline int ldx() const { return ldx_;}
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline Type const& elt2Impl(int i, int j) const  { return p_data()[i*ldx_ + j];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i,j indexes of the element
     **/
    inline Type& elt2Impl(int i, int j) { return p_data()[i*ldx_ + j];}
    /** set a value to this container.
     *  @param v the value to set
     **/
    void setValue(Type const& v)
    {
      for (int i = this->beginRows(); i < this->endRows(); ++i)
        for (int j= this->beginCols(); j < this->endCols(); ++j)
        { this->elt(i, j) = v;}
    }

  protected:
    /** index of the data set */
    int ldx_;
    /** set index of the data. */
    inline void setIdx( int idx) { ldx_ = idx;}
    /** exchange T with this.
     *  @param T the container to move
     **/
    inline void exchange(OrientedCAllocator &T)
    {
      allocator_.exchange(T.allocator_);
      Base::exchange(T);
      std::swap(ldx_, T.ldx_);
    }
    /** @brief Compute the range of the 1D Allocator when we want to
     *  allocate a 2D array with I indexes in the first dimension and J indexes
     *  in the second dimension.
     *  @param I,J range of the rows and columns
     *  @return The range of the 1D allocator
     **/
    static inline Range prod(Range const& I, Range const& J)
    { return Range(J.size()*I.begin()+J.begin(), I.size()*J.size());}
    /** return the increment corresponding to the actual first row an column. */
    inline int shiftInc(int firstRow, int firstCol)
    { return ldx_*firstRow+firstCol; }
    /** set the index corresponding to the actual size of the allocator. */
    inline void setSizedIdx() { ldx_ = this->asDerived().sizeCols();}
    /** manager of the memory */
    Allocator allocator_;
};

/** @ingroup Arrays
 *  @brief  Base class for the general by_col_ structured case.
 **/
template<class Derived, int SizeRows_, int SizeCols_, bool Orient_>
class StructuredCAllocator
    : public OrientedCAllocator<Derived, Orient_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, Orient_> Base;
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J) {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int nbRow, int nbCol)
                               : Base(q, nbRow, nbCol)
    {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T) { return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T) { Base::exchange(T);}

  public:
    /** shift the first indexes of the allocator (for square matrices).
     *  @param firstIdx the index of the first row and column
     **/
    void shift1Impl(int firstIdx)
    { this->asDerived().shift2Impl(firstIdx, firstIdx);}
};

/** @ingroup Arrays
 *  @brief specialization for the point_ case.
 **/
template<class Derived, int SizeCols_>
class StructuredCAllocator<Derived, 1, SizeCols_, Arrays::by_col_>
    : public OrientedCAllocator<Derived, Arrays::by_col_>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, Arrays::by_col_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    using Base::ldx_;
    using Base::p_data;

  protected:
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin()) {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), row_(A.row_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_col_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J), row_(I.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int nbCol)
                               : Base(q, 1, nbCol), row_(0)
    {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T); std::swap(row_, T.row_);}

  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param j index of the column
     **/
    inline Type const& elt1Impl( int j) const { return p_data()[j*ldx_ + row_];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param j index of the columns
     **/
    inline Type& elt1Impl( int j) { return p_data()[j*ldx_ + row_];}
    /** shift the first indexes of the allocator.
     *  @param first the index of the first column and first row */
    inline void shift1Impl(int first)
    { row_ = first; this->asDerived().shift2Impl(first, first);}
    /** resize the allocator.
     *  @param sizeCols the size of the point
     **/
    void resize1Impl(int sizeCols)
    { this->asDerived().resize2Impl(1, sizeCols); row_ = this->beginRows();}
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param J range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& J) const { return Base::row(row_, J);}

  private:
    /** row of the point (needed when this is a reference) */
    int row_;
};

/** @ingroup Arrays
 *  @brief specialization for the point_ case.
 **/
template<class Derived, int SizeCols_>
class StructuredCAllocator<Derived, 1, SizeCols_, Arrays::by_row_>
    : public OrientedCAllocator<Derived, Arrays::by_row_>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, Arrays::by_row_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    using Base::p_data;

  protected:
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin())
                               , p_start_(p_data() + row_*J.size())
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), row_(A.row_)
                               , p_start_(p_data() + row_*A.ldx())
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_row_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J), row_(I.begin())
                               , p_start_(p_data() + row_*A.ldx())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int nbCol)
                               : Base(q, 1, nbCol), row_(0)
                               , p_start_(p_data())
    {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; p_start_ = T.p_start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(p_start_, T.p_start_);
    }
  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param j index of the column
     **/
    inline Type const& elt1Impl( int j) const { return p_start_[j];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param j index of the columns
     **/
    inline Type& elt1Impl( int j) { return p_start_[j];}
    /** shift the first indexes of the allocator.
     *  @param first the index of the first column and first row */
    inline void shift1Impl(int first)
    {
      this->asDerived().shift2Impl(first, first);
      row_ = first;
      p_start_ = p_data() + row_*Base::ldx_;
    }
    /** resize the allocator.
     *  @param sizeCols the size of the point
     **/
    void resize1Impl(int sizeCols)
    { this->asDerived().resize2Impl(1, sizeCols);
      row_ = this->beginRows();
      p_start_ = p_data() + row_*Base::ldx_;
    }
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param J range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& J) const { return Base::row(row_, J);}

  private:
    /** row of the point (needed when this is a reference) */
    int row_;
    /** starting ptr for 1D arrays */
    Type* p_start_;
};

/** @ingroup Arrays
 *  @brief specialization for the vector_ case.
 **/
template<class Derived, int SizeRows_>
class StructuredCAllocator<Derived, SizeRows_, 1, Arrays::by_col_>
    : public OrientedCAllocator<Derived, Arrays::by_col_>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, Arrays::by_col_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    using Base::p_data;

  protected:
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), col_(J.begin())
                               , p_start_(p_data() + col_*I.size())
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), col_(A.col_)
                               , p_start_(p_data() + col_*A.ldx())
    {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_col_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , col_(J.begin())
                               , p_start_(p_data() + col_*A.ldx())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int nbRow, int)
                               : Base(q, nbRow, 1), col_(0)
                               , p_start_(p_data())
                               {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { col_ = T.col_; p_start_ = T.p_start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(col_, T.col_);
      std::swap(p_start_, T.p_start_);
    }

  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type const& elt1Impl( int i) const { return p_start_[i];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type& elt1Impl( int i) { return p_start_[i];}
    /** shift the first indexes of the allocator.
     *  @param first the index of the first row and first column
     **/
    inline void shift1Impl(int first)
    {
      this->asDerived().shift2Impl(first, first);
      col_ = first;
      p_start_ = p_data() + col_*Base::ldx_;
    }
    /** resize the allocator.
     *  @param sizeRow the size of the vector
     **/
    void resize1Impl(int sizeRow)
    {
      this->asDerived().resize2Impl(sizeRow, 1);
      col_ = this->beginCols();
      p_start_ = p_data() + col_*Base::ldx_;
    }
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param I range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& I) const { return Base::col(I, col_);}

  private:
    int col_;
    /** starting ptr for 1D arrays */
    Type* p_start_;
};

/** @ingroup Arrays
 *  @brief specialization for the vector_ case.
 **/
template<class Derived, int SizeRows_>
class StructuredCAllocator<Derived, SizeRows_, 1, Arrays::by_row_>
    : public OrientedCAllocator<Derived, Arrays::by_row_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, Arrays::by_row_> Base;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;
    using Base::ldx_;
    using Base::p_data;

    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), col_(J.begin())
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref), col_(A.col_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_row_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J), col_(J.begin())
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int nbRow, int)
                               : Base(q, nbRow, 1), col_(0)
                               {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { col_ = T.col_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(col_, T.col_);
    }
  public:
    /** @return a constant reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type const& elt1Impl( int i) const { return p_data()[i*ldx_ + col_];}
    /** @return a reference on the element (i,j) of the Allocator.
     *  @param i index of the row
     **/
    inline Type& elt1Impl( int i) { return p_data()[i*ldx_ + col_];}
    /** shift the first indexes of the allocator.
     *  @param firstCol the index of the first column
     **/
    inline void shift1Impl(int firstCol)
    { this->asDerived().shift2Impl(firstCol, firstCol);
      col_ = firstCol;
    }
    /** resize the allocator.
     *  @param sizeRow the size of the vector
     **/
    void resize1Impl(int sizeRow)
    { this->asDerived().resize2Impl(sizeRow, 1); col_ = this->beginCols();}
    /** @return a sub-vector in the specified range of the Allocator.
     *  @param I range of the sub-vector
     **/
    inline SubVector sub1Impl( Range const& I) const { return Base::col(I, col_);}
  private:
    int col_;
};

/** @ingroup Arrays
 *  @brief specialization for the number_ case.
 **/
template<class Derived>
class StructuredCAllocator<Derived, 1, 1, Arrays::by_col_>
    : public OrientedCAllocator<Derived, Arrays::by_col_>
{
  protected:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, Arrays::by_col_> Base;
    using Base::p_data;

    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin()), col_(J.begin())
                               , start_(col_*I.size() + row_)
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref)
                               , row_(A.row_), col_(A.col_)
                               , start_(col_*A.ldx() + row_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_col_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , row_(I.begin()), col_(J.begin())
                               , start_(row_*A.ldx() + col_)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int)
                               : Base(q, 1, 1), row_(0), col_(0)
                               , start_(0)
                               {}
    inline ~StructuredCAllocator() {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; col_ = T.col_; start_ = T.start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(col_, T.col_);
      std::swap(start_, T.start_);
    }
  public:
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt0Impl() const { return p_data()[start_];}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt0Impl() { return p_data()[start_];}
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt1Impl(int) const { return p_data()[start_];}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt1Impl(int) { return p_data()[start_];}

    /** shift the first indexes of the allocator.
     *  @param firstIdx the index of the first row and column
     **/
    inline void shift1Impl(int firstIdx)
    {
      this->asDerived().shift2Impl(firstIdx, firstIdx);
      row_ = firstIdx; col_ = firstIdx; start_ = col_*Base::ldx_ + row_;
    }

  private:
    int row_;
    int col_;
    /** starting idx for number_ arrays */
    int start_;
};

/** @ingroup Arrays
 *  @brief specialization for the number_ case.
 **/
template<class Derived>
class StructuredCAllocator<Derived, 1, 1, Arrays::by_row_>
    : public OrientedCAllocator<Derived, Arrays::by_row_>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef OrientedCAllocator<Derived, Arrays::by_row_> Base;
    using Base::p_data;

  protected:
    /** Default constructor */
    inline StructuredCAllocator( Range const& I, Range const& J)
                               : Base(I, J), row_(I.begin()), col_(J.begin())
                               , start_(row_*J.size() + col_)
    {}
    /** copy constructor */
    inline StructuredCAllocator( StructuredCAllocator const& A, bool ref)
                               : Base(A, ref)
                               , row_(A.row_), col_(A.col_)
                               , start_(row_*A.ldx() + col_) {}
    /** Reference constructor */
    template<class OtherDerived, int OtherSizeRows_, int OtherSizeCols_>
    inline StructuredCAllocator( StructuredCAllocator<OtherDerived, OtherSizeRows_, OtherSizeCols_, Arrays::by_row_> const& A
                               , Range const& I, Range const& J)
                               : Base(A, I, J)
                               , row_(I.begin()), col_(J.begin())
                               , start_(row_*A.ldx() + col_)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline StructuredCAllocator( Type* const& q, int , int)
                               : Base(q, 1, 1), row_(0), col_(0)
                               , start_(0)
                               {}
    /** destructor */
    inline ~StructuredCAllocator() {}
    /** move T to this.
     *  @param T the container to move
     **/
    inline StructuredCAllocator& move(StructuredCAllocator const& T)
    { row_ = T.row_; col_ = T.col_; start_ = T.start_; return *this;}
    /** exchange T with this.
     *  @param T the container to exchange
     **/
    inline void exchange(StructuredCAllocator &T)
    { Base::exchange(T);
      std::swap(row_, T.row_);
      std::swap(col_, T.col_);
      std::swap(start_, T.start_);
    }
  public:
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt0Impl() const { return p_data()[start_];}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt0Impl() { return p_data()[start_];}
    /** @return a constant reference on the element of the Allocator. */
    inline Type const& elt1Impl(int) const { return p_data()[start_];}
    /** @return a reference on the element of the Allocator. */
    inline Type& elt1Impl(int) { return p_data()[start_];}
    /** shift the first indexes of the allocator.
     *  @param firstIdx the index of the first row and column
     **/
    inline void shift1Impl(int firstIdx)
    {
      this->asDerived().shift2Impl(firstIdx, firstIdx);
      row_ = firstIdx; col_ = firstIdx; start_ = col_*Base::ldx_ + row_;
    }
  private:
    int row_;
    int col_;
    /** starting idx for number_ arrays */
    int start_;
};

/** @ingroup Arrays
 *  @brief Allocator for the dense CArray classes.
 *  The size of the Allocator is known in both dimension
 */
template<typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
class CAllocator
      : public StructuredCAllocator<CAllocator<Type_, SizeRows_, SizeCols_, Orient_>, SizeRows_, SizeCols_, Orient_  >
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type, SizeRows_* SizeCols_> Allocator;
    typedef StructuredCAllocator<CAllocator, SizeRows_, SizeCols_, Orient_  > Base;
    using Base::allocator_;

    inline CAllocator(): Base(SizeRows_, SizeCols_) {}
    inline CAllocator( int, int): Base(SizeRows_, SizeCols_) {}
    inline CAllocator( int, int, Type const& v): Base(SizeRows_, SizeCols_) { this->setValue(v);}
    inline CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);} }
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int , int ): Base(q, SizeRows_, SizeCols_) {}
    ~CAllocator() {}
    inline void exchange(CAllocator &T) { Base::exchange(T);}
    inline CAllocator& move(CAllocator const& T)
    {
      if (this == &T) return *this;
      allocator_.move(T);
      Base::move(T);
      Base::setRanges(T.rows(), T.cols());
      Base::setIdx(T.ldx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // check for reference
      if (this->isRef())
      { STKRUNTIME_ERROR_2ARG(CAllocator::shift2Impl, firstRow, firstCol, cannot operate on reference);}
      // set new ranges and translate main pointer
      IContainer2D<SizeRows_, SizeCols_>::shift(firstRow, firstCol);
      allocator_.shiftData(this->shiftInc(firstRow, firstCol));
    }
    inline CAllocator& resize2Impl( int, int) { return *this;}
    inline void realloc(int, int) {}
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The sizes of the columns and of the rows are unknown. The Orientation is
 *  either by rows or by column.
 */
template<typename Type_, bool Orient_>
class CAllocator<Type_, UnknownSize, UnknownSize, Orient_>
     : public StructuredCAllocator<CAllocator<Type_, UnknownSize, UnknownSize, Orient_>, UnknownSize, UnknownSize, Orient_ >
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type, UnknownSize> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, UnknownSize, Orient_ > Base;
    using Base::allocator_;

    /** Default constructor */
    inline CAllocator(): Base(0, 0)  {}
    /** Constructor with specified size.
     *  @param sizeRows, sizeCols size of the rows and columns
     **/
    inline CAllocator( int sizeRows, int sizeCols): Base(sizeRows, sizeCols)
    {}
    /** Constructor with specified size and specified value.
     *  @param sizeRows, sizeCols size of the rows and columns
     *  @param v the initial value
     **/
    inline CAllocator( int sizeRows, int sizeCols, Type const& v)
                     : Base(sizeRows, sizeCols)
    { this->setValue(v);}
    /** Copy or wrapper constructor.
     *  @param A : the array to copy
     *  @param ref : is this a wrapper of A ?
     **/
    inline CAllocator( CAllocator const& A, bool ref = true): Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);}}
    /** Wrapper constructor. This become a reference on (some part of) the Allocator A.
     *  @param A original allocator
     *  @param I,J range of the rows and columns to wrap.
     **/
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int nbRow, int nbCol)
                     : Base(q, nbRow, nbCol)
    {}
    /** Destructor */
    inline ~CAllocator() {}
    /** exchange this with T.
     *  @param T the allocator to exchange
     **/
    inline void exchange(CAllocator &T) { Base::exchange(T);}
    /** move T to this.
     *  @param T the container to move
     **/
    inline CAllocator& move(CAllocator const& T)
    {
      allocator_.move(T.allocator_);
      Base::move(T);
      IContainer2D<UnknownSize, UnknownSize>::setRanges(T.rows(), T.cols());
      this->setIdx(T.ldx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // set new ranges and  translate main pointer
      IContainer2D<UnknownSize, UnknownSize>::shift(firstRow, firstCol);
      allocator_.shiftData(this->shiftInc(firstRow, firstCol));
    }
    CAllocator& resize2Impl( int sizeRows, int sizeCols)
    {
     // check size
     if ((sizeRows <= 0)||(sizeCols<=0))
     {
       // free any allocated memory if this is not a reference
       allocator_.free();
       this->setRanges(sizeRows, sizeCols);
       this->setSizedIdx();
       return *this;
     }
     // allocate
     allocator_.malloc(this->prod(sizeRows, sizeCols));
     this->setRanges(sizeRows, sizeCols);
     this->setSizedIdx();
     return *this;
    }
   /** @brief function for memory reallocation.
    *  The function assume you want to increase or reduce the size without
    *  modification of the bases ranges.
    *  @param sizeRows, sizeCols size of the rows and columns
    **/
   void realloc(int sizeRows, int sizeCols)
   {
     if ((sizeRows == this->sizeRows())&&(sizeCols == this->sizeCols())) return;
     // create a copy the original data set
     CAllocator copy;
     this->exchange(copy);
     try
     {
       // create new container
       resize2Impl(sizeRows, sizeCols);
       shift2Impl(copy.beginRows(), copy.beginCols());
       // copy data
       const int endRow = std::min(copy.endRows(), this->endRows());
       const int endCol = std::min(copy.endCols(), this->endCols());
       for (int j= this->beginCols(); j < endCol; ++j)
         for (int i = this->beginRows(); i< endRow; ++i)
         { this->elt(i, j) = copy.elt(i, j);}
     }
     catch (std::bad_alloc & error)  // if an alloc error occur
     {
       this->exchange(copy); // restore the original container
       STKRUNTIME_ERROR_2ARG(CAllocator::realloc, sizeRows, sizeCols, memory allocation failed);
     }
   }
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The number of rows is known, the number of columns unknown
 **/
template<typename Type_, int SizeRows_, bool Orient_>
class CAllocator<Type_, SizeRows_, UnknownSize, Orient_>
      : public StructuredCAllocator<CAllocator<Type_, SizeRows_, UnknownSize, Orient_>, SizeRows_, UnknownSize, Orient_  >
{
  public:
    typedef Type_ Type;
    typedef StructuredCAllocator<CAllocator, SizeRows_, UnknownSize, Orient_  > Base;
    typedef AllocatorBase<Type, UnknownSize> Allocator;
    using Base::allocator_;

    inline CAllocator(): Base(SizeRows_, 0) {}
    inline CAllocator( int, int sizeCols): Base(SizeRows_, sizeCols) {}
    inline CAllocator( int, int sizeCols, Type const& v): Base(SizeRows_, sizeCols)
    { this->setValue(v);}
    inline CAllocator( CAllocator const& A, bool ref = true)
                     : Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);}}
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int , int nbCol)
                     : Base(q, SizeRows_, nbCol)
    {}
    inline ~CAllocator() {}
    inline void exchange(CAllocator &T) { Base::exchange(T);}
    inline CAllocator& move(CAllocator const& T)
    {
      if (this == &T) return *this;
      allocator_.move(T.allocator_);
      Base::move(T);
      Base::setRanges(T.rows(), T.cols());
      this->setIdx(T.ldx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // set new ranges and  translate main pointer
      IContainer2D<SizeRows_, UnknownSize>::shift(firstRow, firstCol);
      allocator_.shiftData(this->shiftInc(firstRow, firstCol));
    }
    CAllocator& resize2Impl( int, int sizeCols)
    {
      // check size
      if (sizeCols<=0)
      {
        // free any allocated memory if it is not a reference
        allocator_.free();
        this->setRanges(SizeRows_, sizeCols);
        this->setSizedIdx();
        return *this;
      }
      // allocate
      allocator_.malloc(this->prod(SizeRows_, sizeCols));
      this->setRanges(SizeRows_, sizeCols);
      this->setSizedIdx();
      return *this;
    }
    void realloc(int, int sizeCols)
    {     // create a copy the original data set
      CAllocator copy;
      exchange(copy);
      try
      {
        // create new container
        resize2Impl(SizeRows_, sizeCols);
        shift2Impl(copy.beginRows(), copy.beginCols());
        // copy data
        const int endCol = std::min(copy.endCols(), this->endCols());
        for (int j= this->beginCols(); j < endCol; ++j)
          for (int i = this->beginRows(); i < this->endRows(); ++i)
        { this->elt(i, j) = copy.elt(i, j);}

      }
      catch (std::bad_alloc const& error)  // if an alloc error occur
      {
        this->exchange(copy); // restore the original container
        STKRUNTIME_ERROR_2ARG(CAllocator::realloc, SizeRows_, sizeCols, memory allocation failed);
      }
    }
};

/** @brief Specialized Allocator for the dense Arrays classes.
 *  The sizes of the columns is known, the number of rows is unknown
 */
template<typename Type_, bool Orient_, int SizeCols_>
class CAllocator<Type_, UnknownSize, SizeCols_, Orient_>
      : public StructuredCAllocator<CAllocator<Type_, UnknownSize, SizeCols_, Orient_>, UnknownSize, SizeCols_, Orient_  >
{
  public:
    typedef Type_ Type;
    typedef AllocatorBase<Type, UnknownSize> Allocator;
    typedef StructuredCAllocator<CAllocator, UnknownSize, SizeCols_, Orient_  > Base;
    using Base::allocator_;

    inline CAllocator(): Base(0, SizeCols_) {}
    inline CAllocator( int sizeRows, int): Base(sizeRows, SizeCols_) {}
    inline CAllocator( int sizeRows, int, Type const& v)
                     : Base(sizeRows, SizeCols_)
    { this->setValue(v);}
    inline CAllocator( CAllocator const& A, bool ref = true)
                     : Base(A, ref)
    { if (!ref) { allocator_.copy(A.allocator_);} }

    /** wrap other allocator */
    template< int OtherSizeRows_, int OtherSizeCols_>
    inline CAllocator( CAllocator<Type, OtherSizeRows_, OtherSizeCols_, Orient_> const& A
                     , Range const& I, Range const& J)
                     : Base(A, I, J)
    {}
    /** wrapper constructor for 0 based C-Array*/
    inline CAllocator( Type* const& q, int nbRow, int )
                     : Base(q, nbRow, SizeCols_)
    {}
    inline ~CAllocator() {}
    /** exchange T with this without data copy */
    inline void exchange(CAllocator &T) { Base::exchange(T);}
    /** move T to this without data copy */
    inline CAllocator& move(CAllocator const& T)
    {
      if (this == &T) return *this;
      allocator_.move(T.allocator_);
      Base::move(T);
      Base::setRanges(T.rows(), T.cols());
      this->setIdx(T.ldx());
      return *this;
    }
    void shift2Impl(int firstRow, int firstCol)
    {
      // check if there is something to do
      if ((firstRow == this->beginRows())&&(firstCol == this->beginCols())) return;
      // set new ranges and  translate main pointer
      IContainer2D<UnknownSize, SizeCols_>::shift(firstRow, firstCol);
      allocator_.shiftData(this->shiftInc(firstRow, firstCol));
    }
    CAllocator& resize2Impl( int sizeRows, int)
    {
      // check size
      if (sizeRows <= 0)
      {
        // free any allocated memory if it is not a reference
        allocator_.free();
        this->setRanges(sizeRows, SizeCols_);
        this->setSizedIdx();
        return *this;
     }
     // allocate
     allocator_.malloc(this->prod(sizeRows,  SizeCols_));
     this->setRanges(sizeRows,  SizeCols_);
     this->setSizedIdx();
     return *this;
   }
   void realloc(int sizeRows, int)
   {
     // create a copy the original data set
     CAllocator copy;
     exchange(copy);
     try
     {
       // create new container
       resize2Impl(sizeRows, SizeCols_);
       shift2Impl(copy.beginRows(), copy.beginCols());
       // copy data
       const int endRow = std::min(copy.endRows(), this->endRows());
       for (int j= this->beginCols(); j < this->endCols(); ++j)
         for (int i = this->beginRows(); i< endRow; ++i)
         { this->elt(i, j) = copy.elt(i, j);}
     }
     catch (std::bad_alloc & error)  // if an alloc error occur
     {
       this->exchange(copy); // restore the original container
       STKRUNTIME_ERROR_2ARG(CAllocator::realloc, sizeRows, SizeCols_, memory allocation failed);
     }
   }
};

} // namespace STK

#endif /* STK_CALLOCATOR_H */
