/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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

/** @file STK_Law_Geometric.cpp
 *  @brief In this file we implement the Geometric distribution.
 **/

#include "../include/STK_Law_Geometric.h"
#include "Sdk/include/STK_Macros.h"

#ifdef IS_RTKPP_LIB
#include <Rcpp.h>
#endif

namespace STK
{

namespace Law
{

/* constructor */
Geometric::Geometric(Real const& prob) : Base(_T("Geometric")), prob_(prob)
{}
/* destructor */
Geometric::~Geometric() {}

/* @return a geometric random variate . */
Integer Geometric::rand() const
{
#ifdef IS_RTKPP_LIB
  return R::rgeom(prob_);
#else
  return 0;
#endif
}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @return the value of the pdf
 **/
Real Geometric::pdf(Integer const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dgeom(x, prob_, false);
#else
  return 0;
#endif
}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x a binary value
 *  @return the value of the log-pdf
 **/
Real Geometric::lpdf(Integer const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dgeom(x, prob_, true);
#else
  return 0;
#endif
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a Geometric random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real Geometric::cdf(Real const& t) const
{
#ifdef IS_RTKPP_LIB
  return R::pgeom(t, prob_, true, false);
#else
  return 0;
#endif
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
Integer Geometric::icdf(Real const& p) const
{
#ifdef IS_RTKPP_LIB
  return R::qgeom(p, prob_, true, false);
#else
  return 0;
#endif
}

/* @param prob a probability number
 *  @return a Integer random variate.
 **/
Integer Geometric::rand(Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::rgeom(prob);
#else
  return 0;
#endif
}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x a binary value
 *  @param prob a probability number
 *  @return the value of the pdf
 **/
Real Geometric::pdf(Integer x, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::dgeom(x, prob, false);
#else
  return 0;
#endif
}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x a binary value
 *  @param prob a probability number
 *  @return the value of the log-pdf
 **/
Real Geometric::lpdf(Integer x, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::dgeom(x, prob, true);
#else
  return 0;
#endif
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a Geometric random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real Geometric::cdf(Real const& t, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::pgeom(t, prob, true, false);
#else
  return 0;
#endif
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
Integer Geometric::icdf(Real const& p, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::qgeom(p, prob, true, false);
#else
  return 0;
#endif
}


} // namespace Law

} // namespace STK
