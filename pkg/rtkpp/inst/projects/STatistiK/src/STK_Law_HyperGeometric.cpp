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

/** @file STK_Law_HyperGeometric.h
 *  @brief In this file we define the HyperGeometric distribution.
 **/

#include "../include/STK_Law_HyperGeometric.h"
#include "Sdk/include/STK_Macros.h"

#ifdef IS_RTKPP_LIB
#include <Rcpp.h>
#endif

namespace STK
{

namespace Law
{

/* constructor
 *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
 **/
HyperGeometric::HyperGeometric( int nbSuccesses, int nbFailures, int nbDraws)
                              : Base(_T("HyperGeometric"))
                              , nbSuccesses_(nbSuccesses)
                              , nbFailures_(nbFailures)
                              , nbDraws_(nbDraws)
{}
/* destructor */
HyperGeometric::~HyperGeometric() {}

/* @return a random hypergeometric variate. */
Integer HyperGeometric::rand() const
{
#ifdef IS_RTKPP_LIB
  return R::rhyper(nbSuccesses_, nbFailures_, nbDraws_);
#else
  return 0;
#endif
}
/* @brief compute the probability distribution function (density)
 *  Give the value of the pdf at the point x.
 *  @param x an Integer value
 *  @return the value of the pdf
 **/
Real HyperGeometric::pdf(Integer const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dhyper(x, nbSuccesses_, nbFailures_, nbDraws_, false);
#else
  return 0;
#endif
}
/* @brief compute the log probability distribution function
 *  Give the value of the log-pdf at the point x.
 *  @param x an Integer value
 *  @return the value of the log-pdf
 **/
Real HyperGeometric::lpdf(Integer const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dhyper(x, nbSuccesses_, nbFailures_, nbDraws_, true);
#else
  return 0;
#endif
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a HyperGeometric random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real HyperGeometric::cdf(Real const& t) const
{
#ifdef IS_RTKPP_LIB
  return R::phyper(t, nbSuccesses_, nbFailures_, nbDraws_, true, false);
#else
  return 0;
#endif
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
Integer HyperGeometric::icdf(Real const& p) const
{
#ifdef IS_RTKPP_LIB
  return R::qhyper(p, nbSuccesses_, nbFailures_, nbDraws_, true, false);
#else
  return 0;
#endif
}

/* @brief random hypergeometric variate generation.
 *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
 *  @return a Integer random variate.
 **/
Integer HyperGeometric::rand( int nbSuccesses, int nbFailures, int nbDraws)
{
#ifdef IS_RTKPP_LIB
  return R::rhyper(nbSuccesses, nbFailures, nbDraws);
#else
  return 0;
#endif
}
/* @brief compute the probability distribution function.
 *  Give the value of the pdf at the point x.
 *  @param x an Integer value
 *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
 *  @return the value of the pdf
 **/
Real HyperGeometric::pdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws)
{
#ifdef IS_RTKPP_LIB
  return R::dhyper(x, nbSuccesses, nbFailures, nbDraws, false);
#else
  return 0;
#endif
}
/* @brief compute the log probability distribution function.
 *  @param x an Integer value
 *  @param nbSuccesses, nbFailures, nbDraws number of successes, failures, draws
 *  @return the value of the log-pdf
 **/
Real HyperGeometric::lpdf(Integer x, int nbSuccesses, int nbFailures, int nbDraws)
{
#ifdef IS_RTKPP_LIB
  return R::dhyper(x, nbSuccesses, nbFailures, nbDraws, true);
#else
  return 0;
#endif
}
/* @brief compute the cumulative distribution function
 *  Give the probability that a HyperGeometric random variate is less or equal
 *  to t.
 *  @param t a real value
 *  @return the value of the cdf
 **/
Real HyperGeometric::cdf(Real const& t, int nbSuccesses, int nbFailures, int nbDraws)
{
#ifdef IS_RTKPP_LIB
  return R::phyper(t, nbSuccesses, nbFailures, nbDraws, true, false);
#else
  return 0;
#endif
}
/* @brief inverse cumulative distribution function
 *  The quantile is defined as the smallest value @e x such that
 *  <em> F(x) >= p </em>, where @e F is the cumulative distribution function.
 *  @param p a probability number
 **/
Integer HyperGeometric::icdf(Real const& p, int nbSuccesses, int nbFailures, int nbDraws)
{
#ifdef IS_RTKPP_LIB
  return R::qhyper(p, nbSuccesses, nbFailures, nbDraws, true, false);
#else
  return 0;
#endif
}

} // namespace Law

} // namespace STK
