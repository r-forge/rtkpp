/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2015  Serge Iovleff

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
 * Project:  stkpp::STatistiK::Law
 * created on: 23 janv. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_NegativeBinomial.h
 *  @brief In this file we define the NegativeBinomial distribution.
 **/

#include "../include/STK_Law_NegativeBinomial.h"

#ifdef IS_RTKPP_LIB
#include <Rcpp.h>
#endif

namespace STK
{
namespace Law
{

#ifdef IS_RTKPP_LIB

/* @return a Integer random variate . */
inline Integer NegativeBinomial::rand() const
{ return R::rnbinom(size_, prob_);}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @return the value of the pdf
 **/
inline Real NegativeBinomial::pdf(Integer const& x) const
{ return R::dnbinom((double)x, size_, prob_, false);}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x a binary value
 *  @return the value of the log-pdf
 **/
inline Real NegativeBinomial::lpdf(Integer const& x) const
{ return R::dnbinom((double)x, size_, prob_, true);}
/* @brief compute the cumulative distribution function
 *  Give the probability that a NegativeBinomial random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
inline Real NegativeBinomial::cdf(Real const& t) const
{ return R::pnbinom(t, size_, prob_, true, false);}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param prob a probability number
 **/
inline Integer NegativeBinomial::icdf(Real const& p) const
{ return (Integer)R::qnbinom(p, size_, prob_, true, false);}

/* @param prob a probability number
 *  @return a Integer random variate.
 **/
inline Integer NegativeBinomial::rand(int size, Real const& prob)
{ return R::rnbinom(size, prob);}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @param prob a probability number
 *  @return the value of the pdf
 **/
inline Real NegativeBinomial::pdf(Integer x, int size, Real const& prob)
{ return R::dnbinom((double)x, size, prob, false);}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x a binary value
 *  @param prob a probability number
 *  @return the value of the log-pdf
 **/
inline Real NegativeBinomial::lpdf(Integer x, int size, Real const& prob)
{ return R::dnbinom((double)x, size, prob, true);}
/* @brief compute the cumulative distribution function
 *  Give the probability that a NegativeBinomial random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
inline Real NegativeBinomial::cdf(Real const& t, int size, Real const& prob)
{ return R::pnbinom(t, size, prob , true, false);}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param prob a probability number
 **/
inline Integer NegativeBinomial::icdf(Real const& p, int size, Real const& prob)
{ return (Integer)R::qnbinom(p, size, prob , true, false);}

#else

/* @return a Integer random variate . */
Integer NegativeBinomial::rand() const
{
  return 0;
}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @return the value of the pdf
 **/
Real NegativeBinomial::pdf(Integer const& x) const
{
  return 0;
}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x a binary value
 *  @return the value of the log-pdf
 **/
Real NegativeBinomial::lpdf(Integer const& x) const
{
  return 0;
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a NegativeBinomial random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real NegativeBinomial::cdf(Real const& t) const
{
  return 0;
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param prob a probability number
 **/
Integer NegativeBinomial::icdf(Real const& p) const
{
  return 0;
}

/* @param prob a probability number
 *  @return a Integer random variate.
 **/
Integer NegativeBinomial::rand(int size, Real const& prob)
{
  return 0;
}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @param prob a probability number
 *  @return the value of the pdf
 **/
Real NegativeBinomial::pdf(Integer x, int size, Real const& prob)
{
  return 0;
}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x a binary value
 *  @param prob a probability number
 *  @return the value of the log-pdf
 **/
Real NegativeBinomial::lpdf(Integer x, int size, Real const& prob)
{
  return 0;
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a NegativeBinomial random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real NegativeBinomial::cdf(Real const& t, int size, Real const& prob)
{
  return 0;
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param prob a probability number
 **/
Integer NegativeBinomial::icdf(Real const& p, int size, Real const& prob)
{
  return 0;
}

#endif

} // namespace Law

} // namespace STK


