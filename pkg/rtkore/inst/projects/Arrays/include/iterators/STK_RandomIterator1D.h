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
 * created on: 10 mars 2017
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_RandomIterator1D.h
 *  @brief In this file we define and implement the RandomIterator1D and
 *  ConstRandomIterator1D classes.
 **/

#ifndef STK_RANDOMITERATOR1D_H
#define STK_RANDOMITERATOR1D_H

#include "STK_IteratorBase.h"

namespace STK
{
// forward declaration
template<class Derived> struct RandomIterator1D;
template<class Derived> struct ConstRandomIterator1D;

namespace hidden
{
  /** @ingroup hidden
   *  @brief Specialization for the RandomIterator1D iterator class
   **/
  template<class Derived>
  struct IteratorTraits<RandomIterator1D<Derived> >
  {
    typedef std::random_access_iterator_tag iterator_category;
    typedef typename hidden::Traits<Derived>::Type value_type;
    typedef int difference_type; // using only position
    typedef value_type* pointer;
    typedef value_type& reference;
};

  /** @ingroup hidden
   *  @brief Specialization for the RandomIterator1D iterator class
   **/
  template<class Derived>
  struct IteratorTraits<ConstRandomIterator1D<Derived> >
  {
    typedef std::random_access_iterator_tag iterator_category;
    typedef typename hidden::Traits<Derived>::Type value_type;
    typedef int difference_type; // using only position
    typedef value_type const* pointer;
    typedef value_type const& reference;
};

} // namespace hidden

/** @ingroup Arrays
 *  @brief RandomIterator1D allows to loop over the elements of containers derived
 *  from the interface base class STK::ITContainer1D
 **/
template<class Derived>
struct RandomIterator1D: public IteratorBase< RandomIterator1D<Derived> >
{
    typedef  IteratorBase< RandomIterator1D<Derived> > Base;

    typedef typename Base::iterator_category iterator_category;
    typedef typename Base::value_type value_type;
    typedef typename Base::reference reference;
    typedef typename Base::pointer pointer;
    typedef typename Base::difference_type difference_type;

    using Base::pos_;

    // creating
    RandomIterator1D( Derived& array, int pos): Base(pos), array_(array) {}
    RandomIterator1D( RandomIterator1D const& it): Base(it), array_(it.array_) {}
    ~RandomIterator1D() {}
    RandomIterator1D& operator=(RandomIterator1D const& it)
    { array_ = it.array_; pos_= it.pos_; return *this;}

    // getting
    reference operator*()         { return array_[pos_]; }
    pointer operator->()          { return &(array_[pos_]); }
    reference operator[](int pos) { return array_[pos]; }

    // misc
    friend void swap(RandomIterator1D& lhs, RandomIterator1D& rhs)
    {
      Base::swap(lhs, rhs);
      std::swap(lhs.array_, rhs.array_);
    }

  private:
    Derived& array_;
};
template<class Derived>
struct ConstRandomIterator1D: public IteratorBase< ConstRandomIterator1D<Derived> >
{
    typedef  IteratorBase< ConstRandomIterator1D<Derived> > Base;

    typedef typename Base::iterator_category iterator_category;
    typedef typename Base::value_type value_type;
    typedef typename Base::reference reference;
    typedef typename Base::pointer pointer;
    typedef typename Base::difference_type difference_type;

    using Base::pos_;

    // creating
    ConstRandomIterator1D( Derived const& array, int pos)
                         : Base(pos), array_(array) {}
    ConstRandomIterator1D( ConstRandomIterator1D const& it)
                         : Base(it), array_(it.array_) {}
    ~ConstRandomIterator1D() {}
    ConstRandomIterator1D& operator=(ConstRandomIterator1D const& it)
    { array_ = it.array_; pos_= it.pos_; return *this;}
    // getting
    reference operator*() const       { return array_[pos_]; }
    pointer operator->()  const       { return &(array_[pos_]); }
    reference operator[](int pos) const { return array_[pos]; }

    // misc
    friend void swap(ConstRandomIterator1D& lhs, ConstRandomIterator1D& rhs)
    {
      std::swap(lhs.array_, rhs.array_);
      Base::swap(lhs, rhs);
    }

  private:
    Derived const& array_;
};

} // namespace STK

#endif /* STK_RANDOMITERATOR1D_H */
