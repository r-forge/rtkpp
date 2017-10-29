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
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IArray1D.h
 *  @brief Interface base class for the Array1D, this is an internal header file,
 *  included by other Containers library headers.
 *
 *  You should not attempt to use it directly but rather used one of the
 *  derived class like Array1D, except if you want to create your own
 *  Container Class.
 **/

#ifndef STK_IARRAY1D_H
#define STK_IARRAY1D_H

#include "STK_ExprBase.h"
#include "allocators/STK_AllocatorBase.h"
#include "STK_ITContainer1D.h"

namespace STK
{

/** @ingroup Arrays
 *  @brief template one dimensional Array.
 * 
 * An IArray1D is a template one column container implementing the interface
 * base class ITContainer1D.
 **/
template<class Derived >
class IArray1D: public ITContainer1D<Derived>
              , protected AllocatorBase< typename hidden::Traits<Derived>::Type
                                       , hidden::Traits<Derived>::size_
                                       >
{
  public:
    enum
    {
      size_ = hidden::Traits<Derived>::size_
    };
    typedef typename hidden::Traits<Derived>::Type Type;
    typedef typename hidden::Traits<Derived>::RowRange RowRange;
    typedef typename hidden::Traits<Derived>::ColRange ColRange;

  protected:
    typedef AllocatorBase<Type, size_> Allocator;
    typedef ITContainer1D< Derived > Base;

    /** Default constructor. */
    IArray1D();
    /** constructor with a specified Range.
      *  @param I range of the container
     **/
    IArray1D( Range const& I);
    /** Misc constructor with first and last, initialization with a constant.
     *  @param I range of the container
     *  @param v initial value of the container
     **/
    IArray1D( Range const& I, Type const& v);
    /** Copy constructor
     *  @param T the container to copy
     *  @param ref true if T is wrapped
     **/
    IArray1D( const IArray1D &T, bool ref =false);
    /** Copy constructor
     *  @param T the container to copy
     **/
    template<class OtherDerived>
    IArray1D( ExprBase<OtherDerived> const& T);
    /** constructor by reference, ref_=1.
     *  @param T, I the container and the range of data to wrap
     **/
    IArray1D( IArray1D const& T, Range const& I);
    /** constructor by reference, ref_=1.
     *  @param T,I the container and the range of data to wrap
     **/
    template<class OtherDerived>
    IArray1D( IArray1D<OtherDerived> const& T, Range const& I);

    /** destructor: allocated memory is liberated by AllocatorBase base class.*/
    ~IArray1D() {}

  public:
    /** @return @c true if *this is reference container, @c false otherwise */
    inline bool isRef() const { return Allocator::isRef();}
    /** Modify the state of the container: this become a reference (if ref is
     * @c true) or the owner of the data (if ref is @c false).
     * @note To use with care in order to avoid memory leak
     *  @param ref : has top be false if this own its own data
     **/
    void setRef(bool ref) const { Allocator::setRef(ref);}
    /** @return a pointer on the constant data set*/
    inline Type* const& p_data() const { return Allocator::p_data_;}

    /**  @return the range of the rows of the container */
    inline RowRange const& rows() const  { return this->range();}
     /** @return the index of the first element */
    inline int beginRows() const { return this->begin();}
    /**  @return the ending index of the elements */
    inline int endRows() const { return this->end();}
    /**  @return the size of the container */
    inline int sizeRows() const  { return this->size();}

    /** @return the Horizontal range (1 column) */
    inline ColRange cols() const { return ColRange(1);}
    /** @return the index of the first column */
    inline int beginCols() const { return baseIdx;}
    /**  @return the index of the ending column */
    inline int endCols() const  { return baseIdx+1;}
    /** @return the number of columns */
    inline int sizeCols() const  { return 1;};

    /**  @return the index of the last element */
    inline int lastIdxRows() const  { return this->lastIdx();}
    /**  @return the index of the last element */
    inline int lastIdxCols() const  { return baseIdx;}

    /** @return the maximum possible number of elements without reallocation*/
    int capacity() const { return isRef() ? 0 : this->sizeData();}

    /** access to one element.
     *  @param pos index of the element
     **/
    inline Type& elt1Impl(int pos) { return this->data(pos);}
    /** access to one element const.
     *  @param pos index of the const element
     **/
    inline Type const& elt1Impl(int pos) const { return this->data(pos);}
    /** New beginning index for the object.
     *  @param beg the index of the first column to set
     **/
    void shiftImpl(int beg = baseIdx);
    /**  Resize the container.
     * - call @c shift
     * - call @c pushBack if there will be more elements
     * - call @c popBack if three will be less elements
     * @param I the range to set to the Array1D
     **/
    Derived& resizeImpl(Range const& I);
    /** reserve internal memory for at least size elements.
     *  @param size number of elements to reserve
     **/
    void reserve(int size);
    /** Clear the object. Memory is liberated and the
     *  range of the Container is set to 0:-1 or 1:0 (@see baseIdx).
     **/
    void clear();
    /** move T to this.
     *  @note : T is not modified but just set as a reference of the data it was owner.
     *  @param T the container to move to this.
     **/
     void move(Derived const& T);
     /** Add n Elements to the end of the container.
      *  @param n number of elements to add
      **/
     Derived& pushBack( int n=1);
    /** Delete last elts of the container.
     *  @param n number of elts to delete
     **/
     Derived& popBack(int n = 1);
    /** Delete n elements at the pos index to the container.
     *  @param pos index where to delete elements
     *  @param n number of elements to delete (default 1)
    **/
     Derived& erase(int pos, int n=1);
    /** Insert n elements at the position pos of the container. The bound
     *  end_ should be modified at the very end of the insertion as pos
     *  can be a reference to it.
     *  @param pos,n index where to insert the @c n elements (default is 1)
     **/
     Derived& insertElt( int pos, int n =1);
    /** STL compatibility: Insert element @c v in the range @c I of the Array.
     *  @param I range of the index where to insert elements
     *  @param v the value to insert
     **/
     Derived& insert( Range const& I, Type const& v);
    /** STL compatibility: push front an element.
     *  @param v value to append
     **/
     Derived& push_front(Type const& v);
    /** STL compatibility: append an element v.
     *  @param v value to append
     **/
     Derived& push_back(Type const& v);
    /** Swapping the pos1 elt and the pos2 elt.
     *  @param pos1,pos2 positions of the elements to swap
     **/
    void swap(int pos1, int pos2);
    /** exchange this Container with T.
     *  @param T the Array to exchange with this
     **/
    void exchange(IArray1D &T);
    /** overwrite @c this with @c src.
     *  @note If the size match, @c this is not resized, and in this case,
     *  the method take care of the possibly of overlapping.
     *  @param src the container to assign
     **/
    Derived& assign( IArray1D const& src);

    /** set a value to this container.
     *  @param value the value to set
     **/
    Derived& setValue(Type const& value);
  protected:
    /** @return a writable on the constant data set*/
    inline Type* p_data() { return Allocator::p_data_;}
    /** function for memory allocation and initialization.
     *  This method will free all allocated memory owned by this
     *  container before initialization.
     *  @param I range of the container
     **/
    void initialize(Range const& I);
    /** function for memory allocation and initialization.
     *  The range is not set in this method. If a bad_alloc occur, we set the
     *  range of the container to default before throwing it.
     *  @param I range of the container
     **/
    void allocate(Range const& I);
    /** Method for memory deallocation. Memory is liberated and the
     *  range of the Container is set to begin:begin-1.
     **/
    void freeMem();
};

template<class Derived >
IArray1D<Derived>::IArray1D(): Base()
                             , Allocator(Arrays::evalSizeCapacity(0))
{}

template<class Derived >
IArray1D<Derived>::IArray1D( Range const& I)
                          : Base(I)
                           , Allocator(Arrays::evalRangeCapacity(I))
{}

template<class Derived >
IArray1D<Derived>::IArray1D( Range const& I, Type const& v)
                          : Base(I)
                           , Allocator(Arrays::evalRangeCapacity(I))
{ for(int i=this->begin(); i<this->end(); i++) this->data(i) = v;}

/* Copy constructor
 *  @param T the container to copy
 *  @param ref
 **/
template<class Derived >
IArray1D<Derived>::IArray1D( const IArray1D &T, bool ref)
                          : Base(T)
                           , Allocator(T, ref)
{
  if (!ref)
  { allocate(T.range());
    this->memcpy(this->begin(), T, this->range());
  }
}

/* Copy constructor
 *  @param T the container to copy
 **/
template<class Derived >
template<class OtherDerived>
IArray1D<Derived>::IArray1D( ExprBase<OtherDerived> const& T)
                          : Base(T.range())
                           , Allocator(Arrays::evalRangeCapacity(T.range()))
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(OtherDerived);
  for (int i=this->begin(); i<this->end(); i++) this->elt(i)= T.elt(i);
}

template<class Derived >
IArray1D<Derived>::IArray1D( IArray1D const& T, Range const& I)
                          : Base(I)
                           , Allocator(T, true)
{}
/* constructor by reference, ref_=1.
 *  @param T,I the container and the range of data to wrap
 **/
template<class Derived >
template<class OtherDerived>
IArray1D<Derived>::IArray1D( IArray1D<OtherDerived> const& T, Range const& I): Base(T, I) {}

template<class Derived >
void IArray1D<Derived>::shiftImpl(int beg)
{
  // compute increment
  int inc = beg - this->begin();
  if (inc == 0) return;
  // is this structure just a pointer?
  if (this->isRef())
  { STKRUNTIME_ERROR_1ARG(IArray1D::shift,beg,cannot operate on references);}
  // translate cols_ and data
  this->incRange(inc);
  this->shiftData(beg);
}

template<class Derived >
Derived& IArray1D<Derived>::resizeImpl(Range const& I)
{
  // check if there is something to do
  if ( this->range() == I) return this->asDerived();
  if (this->isRef())
  { STKRUNTIME_ERROR_1ARG(IArray1D::resize,I,cannot operate on references);}
  // translate
  shiftImpl(I.begin());
  // compute number of elements to delete or add
  const int inc = I.end() - this->end();
  // adjust size of the container
  if (inc > 0) this->pushBack(inc);  // more elements
  else         this->popBack(-inc);  // less elements
  return this->asDerived();
}

template<class Derived >
void IArray1D<Derived>::reserve(int size)
{
  // nothing to do
  if (size < this->capacity()) return;
  // is this structure a ptr ?
  if (this->isRef())
  { STKRUNTIME_ERROR_1ARG(IArray1D::reserve,size,cannot operate on references);}
  Allocator::realloc(Range(this->begin(), size));
}

template<class Derived >
void IArray1D<Derived>::clear()
{
  if (this->isRef()) return;   // Nothing to do for ref
  this->freeMem();  // Free Mem
  this->setRange(); // Set dimension to default
}

template<class Derived >
void IArray1D<Derived>::move(Derived const& T)
{
  if (this->asPtrDerived() == &T) return; // avoid
  if (!this->isRef()) { freeMem();}
  Allocator::move(T);       // move Allocator part
  this->setRange(T.range());// Set ITContainer1D part
}

template<class Derived >
Derived& IArray1D<Derived>::pushBack( int n)
{
#ifdef STK_ARRAYS_VERY_VERBOSE
    stk_cout << _T("Entering IArray1D<Derived>::pushBack(") << n << _T(")\n");
#endif
  // checks
  if (n <= 0) return this->asDerived();
  if (this->isRef())
  { STKRUNTIME_ERROR_1ARG(IArray1D::pushBack,n,cannot operate on references);}
  // If the container is empty : create it
  if (this->empty()) { initialize(Range(this->begin(), n));}
  else               { insertElt(this->end(), n);}
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::popBack(int n)
{
  // checks
  if (n <= 0) return this->asDerived();
  if (this->isRef())
  { STKRUNTIME_ERROR_1ARG(IArray1D::popBack,n,cannot operate on reference);}
#ifdef STK_BOUNDS_CHECK
  if (this->size()<n)
  { STKOUT_OF_RANGE_1ARG(IArray1D::popBack,n,size() < n);}
#endif
  // update range
  this->decLast(n);
  if (this->size() == 0) this->freeMem(); // release mem if there's no more elts
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::erase(int pos, int n)
{
  // checks
  if (n<=0) return this->asDerived();
  if (this->isRef())
  { STKRUNTIME_ERROR_2ARG(IArray1D::erase,pos, n,cannot operate on reference);}
#ifdef STK_BOUNDS_CHECK
  if (this->begin() > pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::erase,pos, n,begin() > pos);}
  if (this->lastIdx() < pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::erase,pos, n,lastIdx() < pos);}
  if (this->lastIdx() < pos+n-1)
  { STKOUT_OF_RANGE_2ARG(IArray1D::erase,pos, n,lastIdx() < pos+n-1);}
#endif
  // translate remaining elts
  this->memmove(pos, Range(pos+n, this->end() - pos -n));
  //const int last = this->lastIdx()-n;
  //for (int k=pos; k<=last; k++) this->data(k) = this->data(k+n);
  // update dimensions
  this->decLast(n);
  if (this->size() <= 0) this->freeMem(); // if empty release mem
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::insertElt( int pos, int n)
{
  // checks
  if (n <= 0) return this->asDerived();
  if (this->isRef())
  { STKRUNTIME_ERROR_2ARG(IArray1D::insertElt,pos,n,cannot operate on references);}
#ifdef STK_BOUNDS_CHECK
  if (this->begin() > pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::insertElt,pos, n,begin() > pos);}
  if (this->end() < pos)
  { STKOUT_OF_RANGE_2ARG(IArray1D::insertElt,pos, n,end() < pos);}
#endif
  // allocate, if necessary, the memory for the elements
  if (this->capacity() < this->size()+n)
  {
    // temporary container
    IArray1D Taux;
    exchange(Taux);
    // compute range of the container after insertion
    Range range(Taux.range());
    range.incLast(n);
    // initialize
    try
    {
      this->allocate(range);
    }
    catch (runtime_error const& error)   // if an error occur
    {
      exchange(Taux); // restore this
      throw error;    // and send again the Exception
    }
    // set range
    this->setRange(Taux.range());
    // copy original elements
    this->memcpy(Taux.begin(), Taux, Range(Taux.begin(), pos - this->begin()) );
    this->memcpy(pos+n, Taux, Range(pos, this->end()-pos) );
  }
  else // enough space -> shift the last elements
  {
    this->memmove(pos+n, Range(pos, this->end() - pos));
  }
  this->incLast(n);
  return this->asDerived();
}

template<class Derived >
Derived&  IArray1D<Derived>::insert( Range const& I, Type const& v)
{
  this->insertElt(I.begin(), I.size());
  for (int i=I.begin(); i<I.end(); i++) this->data(i) = v;
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::push_front(Type const& v)
{
  insert(Range(this->begin(), 1), v);
  return this->asDerived();
}

template<class Derived >
Derived& IArray1D<Derived>::push_back(Type const& v)
{
  this->pushBack();
  this->back() = v;
  return this->asDerived();
}

template<class Derived >
void IArray1D<Derived>::swap(int pos1, int pos2)
{
#ifdef STK_BOUNDS_CHECK
  if (this->begin() > pos1)
  { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,begin()>pos1);}
  if (this->lastIdx() < pos1)
  { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,lastIdx()<pos1);}
  if (this->begin() > pos2)
  { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,begin()>pos2);}
  if (this->lastIdx() < pos2)
  { STKOUT_OF_RANGE_2ARG(IArray1D::swap,pos1,pos2,lastIdx()<pos2);}
#endif
  // swap
  std::swap(this->data(pos1), this->data(pos2));
}

template<class Derived >
void IArray1D<Derived>::exchange(IArray1D &T)
{
  Allocator::exchange(T);
  Base::exchange(T);
}

template<class Derived >
Derived& IArray1D<Derived>::assign( IArray1D const& src)
{
  if (p_data() == src.p_data())
  {
    if (this->range() != src.range())
    { STKRUNTIME_ERROR_NO_ARG(IArray1D::assign,cannot copy a part of an array on itself);}
    return this->asDerived();
  }
  // Resize if necessary.
  if ( this->sizeRows() != src.sizeRows() ) { this->resize(src.rows());}
  this->memcpy(this->begin(), src, src.range());
  return this->asDerived();
}

/* set a value to this container.
 *  @param value the value to set
 **/
template<class Derived >
Derived& IArray1D<Derived>::setValue(Type const& v)
{
  for(int i=this->begin(); i<this->end(); i++) this->data(i) = v;
  return this->asDerived();
}

template<class Derived >
void IArray1D<Derived>::initialize(Range const& I)
{
  this->clear();
  this->setRef(false);
  this->allocate(I);
  this->setRange(I);
}

template<class Derived >
void IArray1D<Derived>::allocate(Range const& I)
{
  // compute the size necessary (can be 0)
  int size = Arrays::evalSizeCapacity(I.size());
  try
  {
    this->malloc(Range(I.begin(), size));
  }
  catch (runtime_error const& error)
  {
    // if an error occur set default range
    this->setRange();
    throw error;
  }
}

template<class Derived >
void IArray1D<Derived>::freeMem()
{
  // Nothing to do for ref
  if (this->isRef()) return;
  this->free();                 // free allocated memory (Allocator)
  this->setRange(Range(this->begin(), -1)); // set range to default (Base)
}

} // namespace STK

#endif // STK_IARRAY1D_H
