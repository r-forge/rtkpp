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

/** @file STK_SparseAllocator.h
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
 *  and compressed sparse column (CSC) format.
 */
template<typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class SparseAllocator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CAllocator.
 */
template< typename Type_, int SizeRows_, int SizeCols_, bool Orient_>
struct Traits< SparseAllocator<Type_, SizeRows_, SizeCols_, Orient_> >
{
  private:
    class Void {};

  public:
    typedef Type_  Type;
    typedef typename RemoveConst<Type_>::Type const& ReturnType;
    typedef typename RemoveConst<Type_>::Type const& ConstReturnType;

    enum
    {
      orient_    = Orient_,
      sizeRows_  = SizeRows_,
      sizeCols_  = SizeCols_,
      sizeProd_  = ProductSizeRowsBySizeCols<SizeRows_, SizeCols_>::prod_,
      storage_   = Arrays::sparse_ // always sparse
    };

    typedef SparseAllocator<Type_, 1, 1, Orient_> Number;
    typedef SparseAllocator<Type_, 1, SizeCols_, Orient_> Row;
    typedef SparseAllocator<Type_, SizeRows_, 1, Orient_> Col;
    typedef SparseAllocator<Type_, 1, UnknownSize, Orient_> SubRow;
    typedef SparseAllocator<Type_, UnknownSize, 1, Orient_> SubCol;

    /** If one of the Size is 1, we have a Vector (a column) or a Point (a row)
     *  (What to do if both are = 1 : Type_ or array (1,1) ?).
     **/
    typedef typename If< (SizeRows_ == 1)||(SizeCols_ == 1)   // one row or one column
                       , typename If<(SizeCols_ == 1)
                                    , typename If<SizeRows_==1, Number, SubCol>::Result
                                    , SubRow>::Result
                       , Void
                       >::Result SubVector;
    /** Type of the base allocator allocating data */
    typedef Array1D<Type_, sizeProd_> Allocator;
    typedef Type_  Type;

    typedef SparseAllocator<Type_,UnknownSize,UnknownSize, Orient_> SubArray;
};

} // namespace hidden


// forward declaration
template<class Derived, bool Orient_> class OrientedSparseAllocator;

template<class Derived>
class BaseSparseAllocator: ITContainer2D<Derived>
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

    /** @return the vector with the pointers on rows or columns */
    VectorIdx const& ptr() const { return ptr_;}
    /** @return the vecotor with the indexes of the element in the rows or columns */
    VectorIdx const& idx() const { return idx_;}
    /** @return the vector with the non-zero values */
    Array1D<Type> const& data() const { return data_;}

  protected:
    /** @return the vector with the pointers on rows or columns */
    VectorIdx& ptr() { return ptr_;}
    /** @return the vector with the indexes of the element in the rows or columns */
    VectorIdx& idx() { return idx_;}
    /** @return the vector with the non-zero values */
    Array1D<Type>& data() { return data_;}

  private:
    VectorIdx ptr_;
    VectorIdx idx_;
    Array1D<Type> data_;
};

template<class Derived>
class OrientedSparseAllocator<Derived, Arrays::by_col_>: ITContainer2D<Derived>
{
  public:
    typedef Array1D<int>  VectorIdx;
};

template<class Derived>
class OrientedSparseAllocator<Derived, Arrays::by_row_>: ITContainer2D<Derived>
{
  public:
    typedef Array1D<int>  VectorIdx;
};
/** @brief
 *
 */
template<typename Type, int SizeRows_, int SizeCols_, bool Orient_>
class SparseAllocator
{
  public:
    typedef Array1D<int>  VectorIdx;
    typedef Array1D<Type>  elt_vector_type;

    // accessors
    int n_rows() const { return m_n_rows;}
    int n_cols() const { return m_n_cols;}

    int capacity() const { return data_.capacity();}

    VectorIdx const& rowIdx() const { return rowPtr_;}
    VectorIdx const& colIdx() const { return colIdx_;}
    elt_vector_type const& data() const { return data_;}

    int nnz() const = 0;

    /**  */
    std::vector<Type> full () const = 0;

  private:
    SparseAllocator();
    ~SparseAllocator();
    SparseAllocator ( int n_rows, int n_cols, int nz_max = 0)
                    : m_n_rows(n_rows), m_n_cols (n_cols)
    {
      rowPtr_.reserve(nz_max);
      data_.reserve(nz_max);
    }

    int m_n_rows, m_n_cols;

    VectorIdx rowPtr_;
    VectorIdx colIdx_;
    Array1D<Type> data_;

};

} // namespace STK

#endif /* STK_SPARSEALLOCATOR_H */
