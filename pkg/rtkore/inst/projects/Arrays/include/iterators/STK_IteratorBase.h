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
 * created on: 28 marsh 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IteratorBase.h
 *  @brief In this file we base class for Iterators.
 **/

#ifndef STK_ITERATORBASE_H
#define STK_ITERATORBASE_H

#include <iterator>
#include <Sdk/include/STK_IRecursiveTemplate.h>

namespace STK
{

/** @ingroup Arrays
 *  @brief IteratorBase is a base class for all iterators on containers
 *  @tparam Derived the container on which iterate
 **/
template<class Derived>
struct IteratorBase: public IRecursiveTemplate<Derived>
{
    typedef typename hidden::IteratorTraits<Derived>::iterator_category iterator_category;
    typedef typename hidden::IteratorTraits<Derived>::value_type value_type;
    typedef typename hidden::IteratorTraits<Derived>::reference reference;
    typedef typename hidden::IteratorTraits<Derived>::pointer pointer;
    typedef typename hidden::IteratorTraits<Derived>::difference_type difference_type;

  protected:
    /** default constructor.
     *  @param pos the position of the iterator on the array
     **/
    IteratorBase( int pos): pos_(pos) {}
    /** copy constructor.
     *  @param it the iterator to copy
     **/
    IteratorBase( IteratorBase const& it):  pos_(it.pos_) {}
    /** destructor */
    ~IteratorBase() {}

  public:
    // moving
    Derived& operator++()         { ++pos_; return this->asDerived(); }
    Derived& operator++(int junk) { ++pos_; return this->asDerived(); }
    Derived& operator--()         { --pos_; return this->asDerived(); }
    Derived& operator--(int)      { --pos_; return this->asDerived(); }
    Derived& operator+=(int n)    { pos_+=n; return this->asDerived(); }
    Derived& operator-=(int n)    { pos_-=n; return this->asDerived(); }
    friend IteratorBase operator+( IteratorBase const& it, int n)
    { IteratorBase r(it); r+=n ; return r; }
    friend IteratorBase operator+(int n, IteratorBase const& it)
    { IteratorBase r(it); r+=n ; return r; }
    friend IteratorBase operator-( IteratorBase const& it, int n)
    { IteratorBase r(it); r-=n ; return r; }
    friend IteratorBase operator-(int n, IteratorBase const& it)
    { IteratorBase r(it); r-=n ; return r; }
    friend difference_type operator-(IteratorBase it1, IteratorBase it2)
    { return it1.pos_ - it2.pos_;}
    // comparing (only the position is compared)
    bool operator==( IteratorBase const& rhs) { return(pos_ ==rhs.pos_); }
    bool operator!=( IteratorBase const& rhs) { return(pos_!=rhs.pos_); }
    friend bool operator<(IteratorBase const& lhs, IteratorBase const& rhs)
    { return lhs.pos_ < rhs.pos_; };
    friend bool operator>(IteratorBase const& lhs, IteratorBase const& rhs)
    { return lhs.pos_ > rhs.pos_; };
    friend bool operator<=(IteratorBase const& lhs, IteratorBase const& rhs)
    { return lhs.pos_ <= rhs.pos_; };
    friend bool operator>=(IteratorBase const& lhs, IteratorBase const& rhs)
    { return lhs.pos_ >= rhs.pos_; };

    // misc
    friend void swap(IteratorBase& lhs, IteratorBase& rhs)
    { std::swap(lhs.pos_, rhs.pos_);}

  protected:
    int pos_;
};


} // namespace STK

#endif /* STK_ITERATORBASE_H */
