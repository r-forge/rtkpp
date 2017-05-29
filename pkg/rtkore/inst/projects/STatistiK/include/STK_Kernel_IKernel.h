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

/** @file STK_Kernel_IKernel.h
 *  @brief In this file we define the Interface base class for computing a Kernels.
 **/


#ifndef STK_KERNEL_IKERNEL_H
#define STK_KERNEL_IKERNEL_H

#include <Arrays/include/STK_CArraySquare.h>
#include <Sdk/include/STK_IRunner.h>

namespace STK
{

namespace Kernel
{
/** @ingroup Kernel
 *  Interface Base class for the kernels classes.
 */
class IKernel: public IRunnerBase
{
  public:
    /** default constructor */
    inline IKernel(): IRunnerBase(), gram_() {}
    /** copy constructor
     *  @param kernel kernel to copy
     **/
    inline IKernel(IKernel const& kernel): IRunnerBase(kernel), gram_(kernel.gram_) {}
    /** destructor */
    inline virtual ~IKernel() {}
    /** @return the gram matrix
     *  @note if the run method is not used, the gram matrix is empty
     **/
    inline CSquareX const& k() const { return gram_;}
    /** @return the gram matrix (bis)
     *  @note if the run method is not used, the gram matrix is empty
     **/
    inline CSquareX const& gram() const { return gram_;}
    /** @return computed value of the kernel for the ith and jth individuals.
     *  @param i,j indexes of the individuals
     **/
    inline Real kcomp(int i, int j) const { return gram_(i,j);}
    /** @return computed kernel distance between the ith and jth individuals.
     *  @param i,j indexes of the individuals
     **/
    inline Real kdist(int i, int j) const
    { return gram_(i,i)+gram_(j,j)-2*gram_(i,j);}
    /** @return kernel distance between the ith and jth individuals.
     *  @param i,j indexes of the individuals
     **/
    inline Real dist(int i, int j) const
    { return comp(i,i)+comp(j,j)-2*comp(i,j);}

    /** virtual method.
     *  @return diagonale value of the kernel for the ith individuals.
     *  @param i index of the individual
     **/
    virtual inline Real diag(int i) const {return comp(i,i);};
    // pure virtual
    /** pure virtual method.
     *  @return value of the kernel for the ith and jth individuals.
     *  @param i,j indexes of the individuals
     **/
    virtual Real comp(int i, int j) const =0;
    /** pure virtual method.
     *  @return the number of samples (the number of rows in the data set)
     **/
    virtual int nbSample() const =0;
    /** pure virtual method.
     *  @return the number of variables (the number of columns in the data set)
     **/
    virtual int nbVariable() const =0;

  protected:
    /** the resulting gram_ matrix */
    CSquareX gram_;
};


} // namespace Kernel

} // namespace STK

#endif /* STK_KERNEL_IKERNEL_H */
