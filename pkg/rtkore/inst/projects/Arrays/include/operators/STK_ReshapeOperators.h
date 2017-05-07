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
 * Project:  stkpp::Arrays
 * created on: 17 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ReshapeOperators.h
 *  @brief In this file we implement the DiagonalOperator and DiagonalizeOperator classes.
 **/

#ifndef STK_RESHAPEOPERATORS_H
#define STK_RESHAPEOPERATORS_H


#include "Sdk/include/STK_StaticAssert.h"
#include "STK_SlicingOperators.h"

namespace STK
{

// forward declaration
template< typename Array> class DiagonalizeOperator;
template< typename Array> class DiagonalOperator;
template< typename Array> class UpperTriangularizeOperator;
template< typename Array> class LowerTriangularizeOperator;

namespace hidden
{
/** @ingroup hidden
 *  @brief Traits class for DiagonalizeOperator operator
 */
template<typename Lhs>
struct Traits< DiagonalizeOperator <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = ( (Lhs::sizeRows_ != UnknownSize) && (Lhs::structure_!= (int)Arrays::point_) )
                 ?  Lhs::sizeRows_ : UnknownSize,
    sizeCols_  = ( (Lhs::sizeCols_ != UnknownSize) && (Lhs::structure_!= (int)Arrays::vector_) )
                 ?  Lhs::sizeCols_ : UnknownSize,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<DiagonalizeOperator < Lhs> > Row;
  typedef ColOperator<DiagonalizeOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
  typedef typename Lhs::ConstReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for DiagonalOperator operator
 */
template<typename Lhs>
struct Traits< DiagonalOperator <Lhs> >
{
  enum
  {
    structure_ = Arrays::diagonal_,
    orient_    = Lhs::orient_,
    sizeRows_  = ((Lhs::sizeRows_ < Lhs::sizeCols_)) ?  Lhs::sizeRows_ : Lhs::sizeCols_,
    sizeCols_  = sizeRows_,
    storage_   = Lhs::storage_
  };
  typedef RowOperator<DiagonalOperator < Lhs> > Row;
  typedef ColOperator<DiagonalOperator < Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
  typedef typename Lhs::ReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for UpperTriangularizeOperator operator
 */
template<typename Lhs>
struct Traits< UpperTriangularizeOperator<Lhs> >
{
  enum
  {
    structure_ = Arrays::upper_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };

  typedef RowOperator< UpperTriangularizeOperator< Lhs> > Row;
  typedef ColOperator< UpperTriangularizeOperator< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
  typedef typename Lhs::ReturnType ConstReturnType;
};

/** @ingroup hidden
 *  @brief Traits class for UpperTriangularizeOperator operator
 */
template<typename Lhs>
struct Traits< LowerTriangularizeOperator<Lhs> >
{
  enum
  {
    structure_ = Arrays::lower_triangular_,
    orient_    = Lhs::orient_,
    sizeRows_  = Lhs::sizeRows_,
    sizeCols_  = Lhs::sizeCols_,
    storage_   = Lhs::storage_
  };

  typedef RowOperator< LowerTriangularizeOperator< Lhs> > Row;
  typedef ColOperator< LowerTriangularizeOperator< Lhs> > Col;
  typedef typename Lhs::Type Type;
  typedef typename Lhs::ReturnType ReturnType;
  typedef typename Lhs::ReturnType ConstReturnType;
};

} // end namespace hidden

/** @ingroup Arrays
 *  @class DiagonalizeOperator
  *
  * @brief Generic expression when a vector/point expression is "diagonalized".
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * asDiagonal operator.
  *
  * This class represents an expression where a DiagonalizeOperator operator is
  * applied to a vector/point/diagonal expression. It is the return type of the
  * asDiagonal operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalizeOperator type explicitly.
  */
template< typename Lhs>
class DiagonalizeOperator: public ExprBase< DiagonalizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< DiagonalizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::ReturnType ReturnType;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::Row Row;
    typedef typename hidden::Traits< DiagonalizeOperator<Lhs> >::Col Col;
    enum
    {
        structure_ = hidden::Traits< DiagonalizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalizeOperator<Lhs> >::storage_,
        // this is safe as we can use DiagonalizeOperator only on 1D container
        size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;
    /** Constructor */
    inline DiagonalizeOperator( Lhs const& lhs)
                             : Base(), lhs_(lhs)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Lhs);
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().range();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().range();}

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j indexes of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (lhs_.elt(i, j));}
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the ith element
     **/
    inline ConstReturnType elt1Impl(int i) const { return (lhs_.elt(i));}
    /** @return a constant reference on the element of the expression */
    inline ConstReturnType elt0Impl() const { return (lhs_.elt());}

  protected:
    Lhs const& lhs_;
};


/** @ingroup Arrays
 *  @class DiagonalOperator
  *
  * @brief Generic expression when we want to get the diagonal of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * DiagonalOperator operator.
  *
  * This class represents an expression where a diagonal operator is applied to
  * an expression. It is the return type of the diagonal operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name DiagonalOperator type explicitly.
  */
template< typename Lhs>
class DiagonalOperator: public ExprBase< DiagonalOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< DiagonalOperator< Lhs> > Base;
    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::ReturnType ReturnType;
    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::Row Row;
    typedef typename hidden::Traits< DiagonalOperator<Lhs> >::Col Col;
    enum
    {
        structure_ = hidden::Traits< DiagonalOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< DiagonalOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< DiagonalOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< DiagonalOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< DiagonalOperator<Lhs> >::storage_,
        size_      = (sizeRows_ != UnknownSize) ? sizeRows_ : sizeCols_
    };
    /** Type of the Range for the rows */
    typedef TRange<size_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<size_> ColRange;

    /** Constructor */
    inline DiagonalOperator( Lhs const& lhs)
                           : Base(), lhs_(lhs)
                           , range_( lhs_.beginRows(), (size_ != UnknownSize) ? size_ : lhs_.sizeRows())
    {
      if (lhs.rows()!=lhs.cols())
        STKRUNTIME_ERROR_NO_ARG(DiagonalOperatorBase,lhs.rows()!=lhs.cols());
    }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return range_;}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return range_;}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (this->asDerived().lhs().elt(i, j));}
    /** @return a constant reference on the ith element of the expression
     *  @param i index of the ith element
     **/
    inline ConstReturnType elt1Impl(int i) const { return (this->asDerived().lhs().elt(i,i));}
    /** @return a constant reference on the element of the expression */
    inline ConstReturnType elt0Impl() const { return (this->asDerived().lhs().elt());}

  protected:
    Lhs const& lhs_;
    TRange<size_> range_;
};

/** @ingroup Arrays
 *  @class UpperTriangularizeOperator
  *
  * @brief Generic expression when we want to get the upper-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * UpperTriangularizeOperator operator.
  *
  * This class represents an expression where an UpperTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * upper triangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name UpperTriangularizeOperator type explicitly.
  */
template< typename Lhs>
class UpperTriangularizeOperator: public ExprBase< UpperTriangularizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< UpperTriangularizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< UpperTriangularizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< UpperTriangularizeOperator<Lhs> >::ReturnType ReturnType;
    typedef typename hidden::Traits< UpperTriangularizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits< UpperTriangularizeOperator<Lhs> >::Row Row;
    typedef typename hidden::Traits< UpperTriangularizeOperator<Lhs> >::Col Col;
    enum
    {
        structure_ = hidden::Traits< UpperTriangularizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< UpperTriangularizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< UpperTriangularizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< UpperTriangularizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< UpperTriangularizeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline UpperTriangularizeOperator( Lhs const& lhs)
                                     : Base(), lhs_(lhs)
    {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (this->asDerived().lhs().elt(i, j));}

  protected:
    Lhs const& lhs_;
};

/** @ingroup Arrays
 *  @class LowerTriangularizeOperator
  *
  * @brief Generic expression when we want to get the lower-part of a
  * two-dimensional expression.
  *
  * @tparam Lhs the type of the expression to which we are applying the
  * LowerTriangularizeOperator operator.
  *
  * This class represents an expression where an LowerTriangularizeOperator
  * operator is applied to an expression. It is the return type of the
  * lower triangularize operation.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name LowerTriangularizeOperator type explicitly.
  */
template< typename Lhs>
class LowerTriangularizeOperator: public ExprBase< LowerTriangularizeOperator< Lhs> >, public TRef<1>
{
  public:
    typedef ExprBase< LowerTriangularizeOperator< Lhs> > Base;
    typedef typename hidden::Traits< LowerTriangularizeOperator<Lhs> >::Type Type;
    typedef typename hidden::Traits< LowerTriangularizeOperator<Lhs> >::ReturnType ReturnType;
    typedef typename hidden::Traits< LowerTriangularizeOperator<Lhs> >::ConstReturnType ConstReturnType;

    typedef typename hidden::Traits< LowerTriangularizeOperator<Lhs> >::Row Row;
    typedef typename hidden::Traits< LowerTriangularizeOperator<Lhs> >::Col Col;
    enum
    {
        structure_ = hidden::Traits< LowerTriangularizeOperator<Lhs> >::structure_,
        orient_    = hidden::Traits< LowerTriangularizeOperator<Lhs> >::orient_,
        sizeRows_  = hidden::Traits< LowerTriangularizeOperator<Lhs> >::sizeRows_,
        sizeCols_  = hidden::Traits< LowerTriangularizeOperator<Lhs> >::sizeCols_,
        storage_   = hidden::Traits< LowerTriangularizeOperator<Lhs> >::storage_
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** Constructor */
    inline LowerTriangularizeOperator( Lhs const& lhs)
                                     : Base(), lhs_(lhs)
    {}
    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /**  @return the range of the rows */
    inline RowRange const& rowsImpl() const { return lhs().rows();}
    /** @return the range of the Columns */
    inline ColRange const& colsImpl() const { return lhs().cols();}

    /** @return a constant reference on the (i,j) element of the expression.
     *  @param i, j index of the row and of the column
     **/
    inline ConstReturnType elt2Impl(int i, int j) const { return (this->asDerived().lhs().elt(i, j));}

  protected:
    Lhs const& lhs_;
};

} // namespace STK

#endif /* STK_RESHAPEOPERATORS_H */
