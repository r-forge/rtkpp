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
 * Project:  stkpp::Stat::Kernel
 * created on: 5 avr. 2015
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Kernel_IKernelBase.h
 *  @brief In this file we define the Interface base class for computing a Kernels.
 **/


#ifndef STK_KERNEL_IKERNELBASE_H
#define STK_KERNEL_IKERNELBASE_H

#include <Arrays/include/STK_CArraySquare.h>

namespace STK
{
namespace Kernel
{
/** @ingroup Kernel
 *  Interface Base class for the kernels classes.
 */
template<class Array>
class IKernelBase : public IRunnerWithData<Array>
{
  protected:
    /** constructor with a constant pointer on the data set
     *  @param p_data a pointer on a data set that will be "kernelized"
     **/
    inline IKernelBase(Array const* p_data): Base(p_data), gram_() {}
    /** constructor with a constant reference on the data set
     *  @param data a reference on a data set that will be "kernelized"
     **/
    inline IKernelBase(Array const& data): Base(data), gram_() {}

  public:
    typedef IRunnerWithData<Array> Base;
    typedef typename Array::Row RowVector;
    using Base::p_data_;

    /** destructor */
    inline virtual ~IKernelBase() {}
    /** @return the gram matrix*/
    inline CSquareX const& k() const { return gram_;}
    /** @return the gram matrix (bis) */
    inline CSquareX const& gram() const { return gram_;}
    /** Utility method.
     *  @return the computed value of the kernel for the
     *  ith individual and jth individual.
     *  @param i,j indexes of the individuals
     **/
    inline Real kcomp(int i, int j) const { return gram_(i,j);}
    /** Utility method.
     *  @return computed kernel distance between the ith individual and
     *  jth individual using the computed kernel values.
     *  @param i,j indexes of the individuals
     **/
    inline Real kdist(int i, int j) const { return gram_(i,i)+gram_(j,j)-2*gram_(i,j);}

    /** compute the gram matrix. Default implementation using the pure virtual
     *  method k */
    virtual bool run();

    // pure virtuals
    /** compute the kernel between an individual and himself
     *  @param ind individual to evaluate using the kernel*/
    virtual Real kdiag(RowVector const& ind) const = 0;
    /** compute the kernel value between two individuals
     *  @param ind1,ind2 two individuals to evaluate using the kernel metric */
    virtual Real kcomp(RowVector const& ind1, RowVector const& ind2) const = 0;
    /** compute the kernel distance between two individuals
     *  @param ind1,ind2 two individuals to compare using the kernel metric */
    Real kdist(RowVector const& ind1, RowVector const& ind2) const;
  protected:
    /** the resulting gram_ matrix */
    CSquareX gram_;
    /** symmetrize the gram_ matrix using the upper part */
    void symmetrize();
};

template<class Array>
Real IKernelBase<Array>::kdist(RowVector const& ind1, RowVector const& ind2) const
{ return kdiag(ind1)+kdiag(ind2)-2.*kcomp(ind1, ind2);}

template<class Array>
void IKernelBase<Array>::symmetrize()
{
  // lower part
  for (int j= gram_.begin(); j < gram_.end(); ++j)
  {
    for (int i= gram_.begin(); i < j; ++i)
    { gram_(j,i) = gram_(i,j);}
  }
}

template<class Array>
bool IKernelBase<Array>::run()
{
  gram_.resize(p_data_->rows());
  for (int j= gram_.begin(); j < gram_.end(); ++j)
  {
    // create a reference on the current row
    RowVector row_j(p_data_->row(j), true);
    for (int i= gram_.begin(); i < j; ++i)
    { gram_(i,j) = kcomp(p_data_->row(i), row_j);}
    gram_(j, j) = kdiag(row_j);
  }
  // symmetrize
  symmetrize();
  return true;
}

} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_IKERNELBASE_H */
