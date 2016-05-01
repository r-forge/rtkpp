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
 * Project:  stkpp::
 * created on: 5 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_Linear.h
 *  @brief In this file we define the class and methods for computing a Linear Kernel.
 **/


#ifndef STK_KERNEL_LINEAR_H
#define STK_KERNEL_LINEAR_H

#include "STK_Kernel_IKernelBase.h"

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 * The Linear Kernel is a kernel of the form
 * \f[
 * k(x,y) = <x,y>.
 * \f]
 */
template<class Array>
class Linear: public IKernelBase<Array>
{
  public:
    typedef IKernelBase<Array> Base;
    typedef typename Array::Row RowVector;
    using Base::p_data_;
    using Base::gram_;
    using Base::symmetrize;
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     **/
    Linear(Array const* p_data): Base(p_data) {}
    /** constructor with a constant pointer on the data set
     *  @param data a reference on a data set that will be "kernelized"
     **/
    Linear(Array const& data): Base(data) {}
    /** destructor */
    virtual ~Linear() {}
    /** compute the kernel value between two individuals
     *  @param ind1,ind2 two individuals to compare using the kernel metric */
    virtual Real kcomp(RowVector const& ind1, RowVector const& ind2) const;
    /** compute the kernel between an individual and himself
     *  @param ind the individual to evaluate using the kernel
     **/
    virtual Real kdiag(RowVector const& ind) const;
};

template<class Array>
Real Linear<Array>::kcomp(RowVector const& ind1, RowVector const& ind2) const
{ return ind1.dot(ind2);}

template<class Array>
Real Linear<Array>::kdiag(RowVector const& ind) const
{ return ind.norm2();}

} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_LINEAR_H */
