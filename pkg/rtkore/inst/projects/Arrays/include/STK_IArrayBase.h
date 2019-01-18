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
 * created on: 13 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IArrayBase.h
 *  @brief In this file we define the interface base class IIArrayBase
 **/

#ifndef STK_IARRAYBASE_H
#define STK_IARRAYBASE_H

#include "STK_ArrayBase.h"


namespace STK
{

/** @ingroup Arrays
 *  @brief base class for template arrays.
 *
 * This class is the base that is inherited by all containers storing
 * values (matrix, vector, point). Expressions are not arrays. Any derived
 * class can be a lhs in an expression.
 *
 * The common API for these objects is contained in this class.
 *
 *  @tparam Derived is the derived type, e.g., a matrix, vector, point type or
 *  an expression.
 **/
template< class Derived>
class IArrayBase: public ArrayBase<Derived>
{
  public:
    typedef ArrayBase<Derived> Base;
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::TypeConst TypeConst;

    enum
    {
      structure_ = hidden::Traits<Derived>::structure_,
      orient_    = hidden::Traits<Derived>::orient_,
      sizeRows_  = hidden::Traits<Derived>::sizeRows_,
      sizeCols_  = hidden::Traits<Derived>::sizeCols_,
      storage_   = hidden::Traits<Derived>::storage_
    };

  protected:
    /** Default constructor. Default values are cols=(1:0) and rows=(1:0). */
    IArrayBase(): Base() {}
    /** destructor */
    ~IArrayBase() {}

  public:

    // overloaded operators
    /** @return a constant reference on the element (i,j) of the 2D container.
     *  @param i,j row and column indexes
     **/
    inline TypeConst operator()(int i, int j) const
    {
#ifdef STK_BOUNDS_CHECK
       if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, beginRows() > i);}
       if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, endRows() <= i);}
       if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, beginCols() > j);}
       if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, endCols() <= j);}
#endif
      return this->elt(i,j);}
    /** @return a reference on the element (i,j) of the 2D container.
     *  @param i, j indexes of the element to get
     **/
    inline Type& operator()(int i, int j)
    {
#ifdef STK_BOUNDS_CHECK
       if (this->beginRows() > i) { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, beginRows() > i);}
       if (this->endRows() <= i)  { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, endRows() <= i);}
       if (this->beginCols() > j) { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, beginCols() > j);}
       if (this->endCols() <= j)  { STKOUT_OF_RANGE_2ARG(IArrayBase::elt, i, j, endCols() <= j);}
#endif
      return this->elt(i,j);
    }
    /** @return the ith element
     *  @param i index of the element to get
     **/
    inline TypeConst operator[](int i) const
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i) { STKOUT_OF_RANGE_1ARG(IArrayBase::elt, i, begin() > i);}
      if (this->asDerived().end() <= i)  { STKOUT_OF_RANGE_1ARG(IArrayBase::elt, i, end() <= i);}
#endif
      return this->elt(i);
    }
    /** @return a reference on the ith element
     *  @param i index of the element to get
     **/
    inline Type& operator[](int i)
    {
      STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Derived);
#ifdef STK_BOUNDS_CHECK
      if (this->asDerived().begin() > i) { STKOUT_OF_RANGE_1ARG(IArrayBase::elt, i, begin() > i);}
      if (this->asDerived().end() <= i)  { STKOUT_OF_RANGE_1ARG(IArrayBase::elt, i, end() <= i);}
#endif
      return this->elt(i);
    }
    /** @return a constant reference on the number */
    inline TypeConst operator()() const
    {
      return this->elt();
    }
    /** @return the number */
    inline Type& operator()()
    {
      return this->elt();
    }
};

} // namespace STK

#undef MAKE_RESHAPE_OPERATOR

#endif /* STK_ARRAYBASE_H_ */
