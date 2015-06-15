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

/** @file STK_Law_Binomial.cpp
 *  @brief In this file we implement the Binomial distribution.
 **/

#include "../include/STK_Law_Binomial.h"

#ifndef IS_RTKPP_LIB
#include "Analysis/include/STK_Funct_raw.h"
#else
#include <Rcpp.h>
#endif


namespace STK
{

namespace Law
{

#ifdef IS_RTKPP_LIB

inline int Binomial::rand() const { return (int)R::rbinom(n_, prob_);}
inline Real Binomial::pdf(int const& x) const
{ return (Real)R::dbinom((double)x, (double)n_, prob_, false);}
inline Real Binomial::lpdf(int const& x) const
{ return (Real)R::dbinom((double)x, (double)n_, prob_, true);}
inline Real Binomial::cdf(Real const& t) const
{ return (Real)R::pbinom(t, (double)n_, prob_, true, false);}
inline int Binomial::icdf(Real const& p) const
{ return (int)R::qbinom(p, (double)n_, prob_, true, false);}

inline int Binomial::rand(int n, Real const& prob)
{ return (int)R::rbinom(n, prob);}
inline Real Binomial::pdf(int x, int n, Real const& prob)
{ return (Real)R::dbinom(x, (double)n, prob, false);}
inline Real Binomial::lpdf(int x, int n, Real const& prob)
{ return (Real)R::dbinom((double)x, (double)n, prob, true);}
inline Real Binomial::cdf(Real const& t, int n, Real const& prob)
{ return (Real)R::pbinom(t, (double)n, prob, true, false);}
inline int Binomial::icdf(Real const& p, int n, Real const& prob)
{ return (int)R::qbinom(p, (double)n, prob, true, false);}

#else

int Binomial::rand() const
{
  return 0;
}

Real Binomial::pdf(int const& x) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_pdf_raw(x, n_, prob_);
}

Real Binomial::lpdf(int const& x) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_lpdf_raw(x, n_, prob_);
}
Real Binomial::cdf(Real const& t) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(t)) return Arithmetic<Real>::NA();
  return 0.;
}

int Binomial::icdf(Real const& p) const
{
  // trivial cases
  if (Arithmetic<Real>::isNA(p)) return Arithmetic<Real>::NA();
  return 0.;
}

int Binomial::rand(int n, Real const& prob)
{
  return 0;
}

Real Binomial::pdf(int x, int n, Real const& prob)
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_pdf_raw(x, n, prob);
}

Real Binomial::lpdf(int x, int n, Real const& prob)
{
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_lpdf_raw(x, n, prob);
}
Real Binomial::cdf(Real const& t, int n, Real const& prob)
{
  return 0.;
}

int Binomial::icdf(Real const& p, int n, Real const& prob)
{
  return 0.;
}

#endif

} // namespace Law

} // namespace STK

