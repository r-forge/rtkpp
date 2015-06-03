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

/** @file STK_Law_UniformDiscrete.cpp
 *  @brief In this file we implement the uniform law.
 **/

#include "../include/STK_Law_UniformDiscrete.h"
#include "../include/STK_Law_Util.h"

namespace STK
{

namespace Law
{

/* Generate a pseudo UniformDiscrete random variate. */
int UniformDiscrete::rand() const
{ return a_ + int(generator.rand(double(b_ - a_)));}
/* Give the value of the pdf at x.
 *  @param x a real value
 **/
Real UniformDiscrete::pdf( int const& x) const
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a_)||(x > b_)) return 0.;
  return 1./n_;
}
/* Give the value of the log-pdf at x.
 *  @param x a real value
 **/
Real UniformDiscrete::lpdf( int const& x) const
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a_)||(x > b_)) return -Arithmetic<Real>::infinity();
  return -std::log(n_);
}
/* The cumulative distribution function is
 * \f[
 *  F(t; a,b)= \frac{t - a}{b-a}
 * \f]
 *  @param t a real value
 **/
Real UniformDiscrete::cdf( Real const& t) const
{
  if (!Arithmetic<Real>::isFinite(t) ) return t;
  if (t <= a_) return 0.;
  if (t >= b_) return 1.;
  return (b_ - (int)t)/n_;
}

/* The inverse cumulative distribution function is
 * \f[
 * F^{-1}(p; \lambda) = p (b-a) + a.
 * \f]
 *  @param p a probability
 **/
int UniformDiscrete::icdf( Real const& p) const
{
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Exponential::icdf,p,invalid argument);

  if (!Arithmetic<Real>::isFinite(p) ) return p;
  if (p == 1.) return b_;
  if (p == 0.) return a_;
  return(int)((1.-p) * a_ + p * b_);
}

/* Generate a pseudo UniformDiscrete random variate with the specified
 *  parameter.
 *  @param scale the scale of the distribution
 **/
int UniformDiscrete::rand( int a, int b)
{ return a + int(generator.rand(double(b - a)));}
/* Give the value of the pdf at x.
 *  @param x a real value
 *  @param scale the scale of the distribution
 **/
Real UniformDiscrete::pdf( Real const& x, int a, int b)
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a)||(x > b)) return 0.;
  return 1./(b-a);
}
/* Give the value of the log-pdf at x.
 *  @param x a real value
 *  @param scale the scale of the distribution
 **/
Real UniformDiscrete::lpdf( Real const& x, int a, int b)
{
  if (!Arithmetic<Real>::isFinite(x) ) return x;
  if ((x < a)||(x > b)) return -Arithmetic<Real>::infinity();
  return -std::log(b-a);
}

Real UniformDiscrete::cdf(const Real& t, int a, int b)
{
  return (b - t)/(b-a);
}


int UniformDiscrete::icdf(const Real& p, int a, int b)
{ return (int)((1.-p) * a + p * b);}

} // namespace Law

} // namespace STK

