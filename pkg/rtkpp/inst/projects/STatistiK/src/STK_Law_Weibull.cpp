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
 * Purpose:  Weibull probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Weibull.cpp
 *  @brief In this file we implement the Weibull probability distribution.
 **/

#include "../include/STK_Law_Weibull.h"
#ifdef IS_RTKPP_LIB
#include <Rcpp.h>
#endif

namespace STK
{

namespace Law
{

#ifdef IS_RTKPP_LIB

/*inline*/ Real Weibull::rand() const
{ return R::rweibull(k_, lambda_);}
/*inline*/ Real Weibull::pdf(Real const& x) const
{ return R::dweibull(x, k_, lambda_, false);}
/*inline*/ Real Weibull::lpdf(Real const& x) const
{ return R::dweibull(x, k_, lambda_, true);}
/*inline*/ Real Weibull::cdf(Real const& t) const
{ return R::pweibull(t, k_, lambda_, true, false);}
/*inline*/ Real Weibull::icdf(Real const& p) const
{ return R::qweibull(p, k_, lambda_, true, false);}

/*inline*/ Real Weibull::rand( Real const& k, Real const& lambda)
{ return R::rweibull(k, lambda);}
/*inline*/ Real Weibull::pdf(Real const& x, Real const& k, Real const& lambda)
{ return R::dweibull(x, k, lambda, false);}
/*inline*/ Real Weibull::lpdf(Real const& x, Real const& k, Real const& lambda)
{ return R::dweibull(x, k, lambda, true);}
/*inline*/ Real Weibull::cdf(Real const& t, Real const& k, Real const& lambda)
{ return R::pweibull(t, k, lambda, true, false);}
/*inline*/ Real Weibull::icdf( Real const& p, Real const& k, Real const& lambda)
{ return R::qweibull(p, k, lambda, true, false);}

#else

/* @return a pseudo Weibull random variate. */
Real Weibull::rand() const
{
  return 0;
}
/* @return the value of the pdf
 *  @param x a positive real value
 **/
Real Weibull::pdf(Real const& x) const
{
  return 0;
}
/* @return the value of the log-pdf
 *  @param x a positive real value
 **/
Real Weibull::lpdf(Real const& x) const
{
  return 0;
}
/*The cumulative distribution function for the Weibull distribution is
 *  \f$ F(x;k,\lambda) = 1- e^{-(x/\lambda)^k}.\f$
 *  @return the cumulative distribution function
 *  @param t a positive real value
 **/
Real Weibull::cdf(Real const& t) const
{
  return 0;
}
/*The quantile (inverse cumulative distribution) function for the Weibull
 * distribution is \f$ Q(p;k,\lambda) = \lambda {(-\ln(1-p))}^{1/k} \f$
 *  @return the inverse cumulative distribution function
 *  @param p a probability number
 **/
Real Weibull::icdf(Real const& p) const
{
  return 0;
}

/* @return a pseudo Weibull random variate with the specified parameters.
 *  @param k, lambda shape and scale (dispersion) parameters
 **/
Real Weibull::rand( Real const& k, Real const& lambda)
{
  return 0;
}
/* @return the value of the pdf
 *  @param x a positive real value
 *  @param k, lambda shape and scale (dispersion) parameters
 **/
Real Weibull::pdf(Real const& x, Real const& k, Real const& lambda)
{
  return 0;
}
/* @return the value of the log-pdf
 *  @param x a positive real value
 *  @param k, lambda shape and scale (dispersion) parameters
 **/
Real Weibull::lpdf(Real const& x, Real const& k, Real const& lambda)
{
  return 0;
}
/* @return the cumulative distribution function
 *  @param t a positive real value
 *  @param k, lambda shape and scale (dispersion) parameters
 **/
Real Weibull::cdf(Real const& t, Real const& k, Real const& lambda)
{
  return 0;
}
/* @brief Compute the inverse cumulative distribution function at p
 *  of the standard log-normal distribution.
 *
 *  @param p a probability number.
 *  @param mu, sigma location and scale of the log-normal law
 *  @return the inverse cumulative distribution function value at p.
 **/
Real Weibull::icdf( Real const& p, Real const& k, Real const& lambda)
{
  return 0;
}

#endif

} // namespace Law

} // namespace STK
