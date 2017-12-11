/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2017  Serge Iovleff, Université Lille 1, Inria

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
 * created on: 10 mars 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_SAllocator.h
 *  @brief In this file we implement the sparse allocators.
 **/

#ifndef STK_SPARSEALLOCATOR_H
#define STK_SPARSEALLOCATOR_H

#include "../STK_Array1D.h"

namespace STK
{
/** @ingroup Arrays
 *  @brief Allocator for sparse Array classes.
 *  The data are stored either in the compressed sparse row (CSR) format
 *  or compressed sparse column (CSC) format.
 */
template<typename Type, int Size_, bool Orient_>
class SAllocator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CAllocator.
 */
template< typename Type_, int Size_, bool Orient_>
struct Traits< SAllocator<Type_, Size_, Orient_> >
{
  private:
    class Void {};

  public:
    typedef Type_  Type;
    typedef typename RemoveConst<Type_>::Type const& ConstReturnType;

    enum
    {
      orient_    = Orient_,
//      sizeRows_  = SizeRows_,
//      sizeCols_  = SizeCols_,
      storage_   = Arrays::sparse_ // always sparse
    };

    typedef SAllocator<Type_, 1, Orient_> Number;
//    typedef SAllocator<Type_, 1, SizeCols_, Orient_> Row;
//    typedef SAllocator<Type_, SizeRows_, 1, Orient_> Col;
//    typedef SAllocator<Type_, 1, UnknownSize, Orient_> SubRow;
//    typedef SAllocator<Type_, UnknownSize, 1, Orient_> SubCol;

    /** If one of the Size is 1, we have a Vector (a column) or a Point (a row)
     *  (What to do if both are = 1 : Type_ or array (1,1) ?).
     **/
//    typedef typename If< (SizeRows_ == 1)||(SizeCols_ == 1)   // one row or one column
//                       , typename If<(SizeCols_ == 1)
//                                    , typename If<SizeRows_==1, Number, SubCol>::Result
//                                    , SubRow>::Result
//                       , Void
//                       >::Result SubVector;

    typedef SAllocator<Type_,UnknownSize, Orient_> SubArray;
};

} // namespace hidden


// forward declaration
template<class Derived, bool Orient_> class OrientedSAllocator;

/** @ingroup Arrays base class for spase allocators */
template<class Derived>
class SAllocatorBase
{
  public:

    typedef typename hidden::Traits<Derived>::Type Type;

    typedef Array1D<int>  VectorIdx;

    /** Type for the position and value of an element in a row (or column) */
    typedef std::pair<int, Type> IndexedValue;
    /** Vector storing the elements in a row (or column) */
    typedef Array1D<IndexedValue>  VectorIdxValue;

    /** constructor with specified dimensions
     * @param dim number of column
     * @param nz number of data
     **/
    SAllocatorBase( int dim, int nz)
                  : ptr_(dim), VectorIdxValue(nz)
    {}
    /** copy constructor */
    SAllocatorBase( SAllocatorBase const& A, bool ref)
                  : ptr_(A.ptr_, ref), VectorIdxValue(A, ref)
                  , idx_(A.idx_, ref), data_(A.data_, ref)
    {}

    /** @return the vector with the pointers on rows or columns */
    VectorIdx const& ptr() const { return ptr_;}
    /** @return the vecotor with the indexes of the element in the rows or columns */
    VectorIdx const& idx() const { return idx_;}
    /** @return the vector with the non-zero values */
    Array1D<Type> const& data() const { return data_;}

    /** This method allows to get the element (p_idx, s_idx)
     *  @param p_idx the index of the row (respectively column)
     *  @param s_idx the index of the column (respectively row)
     *  @return 0 if the element is not stored, the value of the element otherwise
     **/
    Type getElt(int p_idx, int s_idx) const
    {
      for (int t=ptr_[p_idx]; t<ptr_[p_idx+1]; ++t)
      { if (idxValues_[t].first == s_idx) return idxValues_[t].second;}
      return Type(0);
    }
    /** This method allows to overwrite or insert an element to the position (p_idx, s_idx)
     *  @param p_idx index of the row (respectively column)
     *  @param s_idx index of the column (respectively row)
     *  @param value value to set
     **/
    void setElt(int p_idx, int s_idx, Type const& value)
    {
      for (int t=ptr_[p_idx]; t<ptr_[p_idx+1]; ++t)
      {
        // existing value for this entry
        if (idxValues_[t].first == s_idx) { idxValues_[t].second = value; return;}
        // we
        if (idxValues_[t].first > s_idx)
        { idxValues_[t].insert(std::pair(s_idx, value)); return;}
      }
      idxValues_.push_back(std::pair(s_idx, value));
    }

  protected:
    /** @return the vector with the pointers on rows (respectively columns) */
    VectorIdx& ptr() { return ptr_;}
    /** @return the vector with the indexes of the element in the row (respectively columns) */
    VectorIdx& idx() { return idx_;}
    /** @return the vector with the non-zero values */
    Array1D<Type>& data() { return data_;}

  private:
    VectorIdx ptr_;
    VectorIdxValue idxValues_;
    // will be needed for external usage
    VectorIdx idx_;
    Array1D<Type> data_;
};

/** Specialization for the row oriented storage scheme */
template<class Derived>
class OrientedSAllocator<Derived, Arrays::by_col_>: ITContainer2D<Derived>
{
  public:
    typedef ITContainer2D<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;
    typedef typename hidden::Traits<Derived>::SubRow SubRow;
    typedef typename hidden::Traits<Derived>::SubCol SubCol;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    typedef Array1D<int>  VectorIdx;

    /** constructor with specified ranges */
    OrientedSAllocator( Range const& I, Range const& J)
                     : Base(I, J)
    {}
    /** copy constructor */
    OrientedSAllocator( OrientedSAllocator const& A, bool ref)
                     : Base(A)
    { if (!ref) allocator_.assign(A.allocator_);}
    /** Reference constructor */
    template<class OtherDerived>
    inline OrientedSAllocator( OrientedCAllocator<OtherDerived, Arrays::by_row_> const& A
                             , Range const& I, Range const& J)
                            : Base(I, J), ldx_(A.ldx()), allocator_(A.allocator(), true)
    {}

    /** @return the vector with the pointers on rows or columns */
    VectorIdx const& col_ptr() const { return ptr_;}
    /** @return the vector with the indexes of the element in the rows or columns */
    VectorIdx const& idx() const { return idx_;}
    /** @return the vector with the non-zero values */
    Array1D<Type> const& data() const { return data_;}

};

template<class Derived>
class OrientedSAllocator<Derived, Arrays::by_row_>: ITContainer2D<Derived>
{
  public:
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::Row Row;
    typedef typename hidden::Traits<Derived>::Col Col;
    typedef typename hidden::Traits<Derived>::SubRow SubRow;
    typedef typename hidden::Traits<Derived>::SubCol SubCol;
    typedef typename hidden::Traits<Derived>::SubArray SubArray;
    typedef typename hidden::Traits<Derived>::SubVector SubVector;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;
    typedef Array1D<int>  VectorIdx;

};
/** @brief
 *
 */
template<typename Type, int Size_, bool Orient_>
class SAllocator: public OrientedSAllocator<SAllocator<Type, Size_, Orient_>, Orient_>
{
  public:
    ~SAllocator();
    SAllocator ( int n_rows, int n_cols, int nz_max = 0)
               : m_n_rows(n_rows), m_n_cols (n_cols)
    {
    }
};

} // namespace STK

#endif /* STK_SPARSEALLOCATOR_H */
