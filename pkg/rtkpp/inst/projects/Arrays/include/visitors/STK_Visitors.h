/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2012  Serge Iovleff

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
 * created on: 27 sept. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Visitors.h
 *  @brief In this file we define the Visitors classes.
 **/

#ifndef STK_VISITORS_H
#define STK_VISITORS_H

#include "STatistiK/include/STK_Law_Util.h"

#include "STK_VisitorsImpl.h"
#include "STK_VisitorSelector.h"


namespace STK
{

namespace hidden
{

/** @ingroup hidden
  * @brief utility class for getting the result from a visitor acting on a
  * vector or a point.
  */
template <class Visitor, int Structure_> struct GetIdx;
template <class Visitor> struct GetIdx<Visitor, Arrays::vector_>
{ static int idx(Visitor const& visitor) { return visitor.row_;};};

template <class Visitor> struct GetIdx<Visitor, Arrays::point_>
{ static int idx(Visitor const& visitor) { return visitor.col_;};};

template <class Visitor> struct GetIdx<Visitor, Arrays::diagonal_>
{ static int idx(Visitor const& visitor) { return visitor.col_;};};


/** @ingroup hidden
  * @brief Base class to implement min, max, sum,... visitors for 2D containers.
  */
template <typename Type>
struct EltVisitor2DBase
{
  typedef Type return_type;
  EltVisitor2DBase() : row_(baseIdx), col_(baseIdx), res_(Arithmetic<Type>::NA()) {};
  int row_, col_;
  Type res_;
  return_type result() const { return res_;}
};

/** @ingroup hidden
  * @brief Visitor computing the min coefficient with its value and coordinates
  *
  * @sa STK::ExprBase::minElt(int, int), STK::ExprBase::minElt(int)
  */
template <typename Type>
struct MinEltVisitor : EltVisitor2DBase<Type>
{
  MinEltVisitor() : EltVisitor2DBase<Type>()
  { this->res_ = Arithmetic<Type>::max(); }
  void operator()( Type const& value, int i, int j)
  {
    if (value < this->res_)
    { this->res_ = value; this->row_ = i; this->col_ = j;}
  }
};

/** @ingroup hidden
  * @brief Visitor computing safely the minimal coefficient with its value and
  * indexes.
  *
  * @sa STK::ExprBase::minSafeElt(int, int), STK::ExprBase::minSafeElt(int)
  */
template <typename Type>
struct MinEltSafeVisitor : EltVisitor2DBase<Type>
{
  MinEltSafeVisitor() : EltVisitor2DBase<Type>()
  { this->res_ = Arithmetic<Type>::max(); }
  void operator()( Type const& value, int i, int j)
  {
    if (Arithmetic<Type>::isFinite(value) && (value < this->res_))
    { this->res_ = value; this->row_ = i; this->col_ = j;}
  }
};
/** @ingroup hidden
 *  @brief Visitor computing the maximal coefficient of the Array
 *
 * @sa ExprBase::maxElt(int, int)
 */
template <typename Type>
struct MaxEltVisitor : EltVisitor2DBase<Type>
{
  MaxEltVisitor(): EltVisitor2DBase<Type>()
  { this->res_ = -Arithmetic<Type>::max(); }
  void operator() ( Type const& value, int i, int j)
  {
    if (value > this->res_)
    { this->res_ = value; this->row_ = i; this->col_ = j;}
  }
};
/** @ingroup hidden
 *  @brief Visitor computing safely the maximal coefficient of the Array
 *
 * @sa ExprBase::maxElt(int, int)
 */
template <typename Type>
struct MaxEltSafeVisitor : EltVisitor2DBase<Type>
{
  MaxEltSafeVisitor(): EltVisitor2DBase<Type>()
  { this->res_ = -Arithmetic<Type>::max(); }
  void operator() ( Type const& value, int i, int j)
  {
    if (Arithmetic<Type>::isFinite(value))
      if(value > this->res_)
      { this->res_ = value; this->row_ = i; this->col_ = j;}
  }
};
/** @ingroup hidden
 *  @brief Visitor computing the min of all the coefficients of the Array
 *
 * @sa ExprBase::minElt(int, int)
 */
template <typename Type>
struct MinVisitor
{
  typedef Type return_type;
  Type res_;
  MinVisitor(): res_(Arithmetic<Type>::max()) {}
  void operator() ( Type const& value, int i, int j)
  { res_ = std::min(res_,value);}
  void operator() ( Type const& value, int i)
  { res_ = std::min(res_,value);}
  return_type result() const { return res_;}
};
/** @ingroup hidden
 *  @brief Visitor computing safely the min of all the coefficients of the Array
 *
 * @sa ExprBase::minElt(int, int)
 */
template <typename Type>
struct MinSafeVisitor
{
  typedef Type return_type;
  Type res_;
  MinSafeVisitor(): res_(Arithmetic<Type>::max()) {}
  void operator() ( Type const& value, int i, int j)
  {  if (Arithmetic<Type>::isFinite(value))
    {res_ = std::min(res_,value);}
  }
  void operator() ( Type const& value, int i)
  { if (Arithmetic<Type>::isFinite(value))
    { res_ = std::min(res_,value);}
  }
  return_type result() const { return res_;}
};
/** @ingroup hidden
 *  @brief Visitor computing the max of all the coefficients of the Array
 *
 * @sa ExprBase::maxElt(int, int), ExprBase::min()
 */
template <typename Type>
struct MaxVisitor
{
  typedef Type return_type;
  Type res_;
  MaxVisitor(): res_(-Arithmetic<Type>::max()) {}
  void operator() ( Type const& value, int i, int j)
  { res_ = std::max(res_,value);}
  void operator() ( Type const& value, int i)
  { res_ = std::max(res_,value);}
  return_type result() const { return res_;}
};
/** @ingroup hidden
 *  @brief Visitor computing safely the max of all the coefficients of the Array
 *
 * @sa ExprBase::max(), ExprBase::min()
 */
template <typename Type>
struct MaxSafeVisitor
{
  typedef Type return_type;
  Type res_;
  MaxSafeVisitor(): res_(-Arithmetic<Type>::max()) {}
  void operator() ( Type const& value, int i, int j)
  { if (Arithmetic<Type>::isFinite(value))
    { res_ = std::max(res_,value);}
  }
  void operator() ( Type const& value, int i)
  { if (Arithmetic<Type>::isFinite(value))
    { res_ = std::max(res_,value);}
  }
  return_type result() const { return res_;}
};
/** @ingroup hidden
 *  @brief Visitor computing the sum of all the coefficients of the Array
 *
 * @sa ExprBase::sum()
 */
template <typename Type>
struct SumVisitor
{
  typedef Type return_type;
  Type res_;
  SumVisitor(): res_(Type(0)) {}
  void operator() ( Type const& value, int i, int j)
  { res_ += value;}
  void operator() ( Type const& value, int i)
  { res_ += value;}
  return_type result() const { return res_;}
};
/** @ingroup hidden
 *  @brief Visitor computing the mean of all the coefficients of the Array
 *
 * @sa ExprBase::mean()
 */
template <typename Type>
struct MeanVisitor
{
  typedef Type return_type;
  Type res_;
  int nb_;
  MeanVisitor(): res_(Type(0)), nb_(0) {}
  void operator() ( Type const& value, int i, int j)
  { res_ += value; nb_++;}
  void operator() ( Type const& value, int i)
  { res_ += value; nb_++;}
  return_type result() const { return nb_ == 0  ? Arithmetic<Type>::NA() : res_/nb_;}
};
/** @ingroup hidden
 *  @brief Visitor computing safely the mean of all the coefficients of the Array
 *
 * @sa ExprBase::mean()
 */
template <typename Type>
struct MeanSafeVisitor
{
  typedef Type return_type;
  Type res_;
  int nb_;
  MeanSafeVisitor(): res_(Type(0)), nb_(0) {}
  void operator() ( Type const& value, int i, int j)
  {
    if (Arithmetic<Type>::isFinite(value))
    {  res_ += value; nb_++;}
  }
  void operator() ( Type const& value, int i)
  {
    if (Arithmetic<Type>::isFinite(value))
    {  res_ += value; nb_++;}
  }
  return_type result() const { return nb_ == 0  ? Arithmetic<Type>::NA() : res_/nb_;}
};
/** @ingroup hidden
 *  @brief Visitor counting the number of not-zero element in an array
 *  This visitor can be used in conjunction with the comparison operators
 *  in order to get the number of element matching a condition. For example:
 *  @code
 *    // get in c the number of element of A equal to 2
 *    int c = (A == 2).count()
 *  @endcode
 *
 *  @sa ExprBase::count()
 */
template <typename Type>
struct CountVisitor
{
  typedef int return_type;
  int res_;
  CountVisitor(): res_(0) {}
  void operator() ( Type const& value, int i, int j)
  { if (value) ++res_;}
  void operator() ( Type const& value, int i)
  { if (value) ++res_;}
  return_type result() const { return res_;}
};

/** @ingroup hidden
 *  @brief Visitor checking if all the elements of an array are different
 *  from zero. For example:
 *  @code
 *    // check if all the elements of A equal to 2
 *    bool f = (A == 2).all()
 *  @endcode
 *
 *  @sa ExprBase::all()
 */
template <typename Type>
struct AllVisitor
{
  typedef bool return_type;
  bool res_;
  AllVisitor(): res_(true) {}
  void operator() ( Type const& value, int i, int j)
  { res_ &= (value);}
  void operator() ( Type const& value, int i)
  { res_ &= (value);}
  return_type result() const { return res_;}
};
/** @ingroup hidden
 *  @brief Visitor checking if at least, one element of an array is different
 *  from zero. For example:
 *  @code
 *    // check if one of the elements of A is equal to 2
 *    bool f = (A == 2).any()
 *  @endcode
 *
 *  @sa ExprBase::all()
 */
template <typename Type>
struct AnyVisitor
{
  typedef bool return_type;
  bool res_;
  AnyVisitor(): res_(false) {}
  void operator() ( Type const& value, int i, int j)
  { res_ |= (value);}
  void operator() ( Type const& value, int i)
  { res_ |= (value);}
  return_type result() const { return res_;}
};

/** @ingroup hidden
 *  @brief Applier setting uniform random numbers
 *
 * @sa STK::ExprBase::randGauss(), STK::ExprBase::rand()
 */
template <typename Type>
struct RandUnifApplier
{
  void operator() ( Type& value)
  { value = Type(Law::generator.randUnif());}
};

/** @ingroup hidden
 *  @brief Visitor setting Gaussian random variables
 *  @sa STK::ExprBase::RandUnif()
 */
template <typename Type>
struct RandGaussApplier
{
  void operator() ( Type& value)
  { value = Type(Law::generator.randGauss());}
};

/** @ingroup hidden
 *  @brief Visitor putting a choosen random variable
 */
template <typename Type>
struct RandApplier
{
  RandApplier( Law::IUnivLaw<Type> const& law):law_(law){}
  void operator() ( Type& value)
  { value = law_.rand();}
  Law::IUnivLaw<Type> const& law_;
};

/** @ingroup hidden
  * @brief visitor putting a constant value
  */
template <typename Type>
struct ValueApplier
{
  Type value_;
  ValueApplier(Type const& value) : value_(value) {};
  void operator() ( Type& value)
  { value = value_;}
};

} //namespace hidden

} // namespace STK

#endif /* STK_VISITORS_H */
