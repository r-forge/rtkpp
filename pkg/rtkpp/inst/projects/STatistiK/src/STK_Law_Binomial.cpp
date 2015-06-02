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

/** @file STK_Law_Binomial.cpp
 *  @brief In this file we implement the Binomial distribution.
 **/

#include "../include/STK_Law_Binomial.h"
#ifdef IS_RTKPP_LIB
#include <Rcpp.h>
#else
#include "Analysis/include/STK_Funct_raw.h"
#endif

namespace STK
{

namespace Law
{
Binomial::Binomial( int n, Real const& prob)
                  : Base(_T("Binomial")), n_(n), prob_(prob)
{
  if (prob<0) STKDOMAIN_ERROR_2ARG(Binomial::Binomial,prob,n,prob must be >= 0);
  if (prob>1) STKDOMAIN_ERROR_2ARG(Binomial::Binomial,prob,n,prob must be <= 1);
  if (n<0) STKDOMAIN_ERROR_2ARG(Binomial::Binomial,prob,n,n must be >= 0);
}
Binomial::~Binomial() {}

int Binomial::rand() const
{
#ifdef IS_RTKPP_LIB
  return (int)R::rbinom(n_, prob_);
#else
  return 0;
#endif
}

Real Binomial::pdf(int const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dbinom((double)x, (double)n_, prob_, false);
#else
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_pdf_raw(x, n_, prob_);
#endif
}

Real Binomial::lpdf(int const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dbinom((double)x, (double)n_, prob_, true);
#else
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_lpdf_raw(x, n_, prob_);
#endif
}
Real Binomial::cdf(Real const& t) const
{
#ifdef IS_RTKPP_LIB
  return R::pbinom(t, (double)n_, prob_, true, false);
#else
  // trivial cases
  if (Arithmetic<Real>::isNA(t)) return Arithmetic<Real>::NA();
  return 0.;
#endif
}

int Binomial::icdf(Real const& p) const
{
#ifdef IS_RTKPP_LIB
  return (int)R::qbinom(p, (double)n_, prob_, true, false);
#else
  // trivial cases
  if (Arithmetic<Real>::isNA(p)) return Arithmetic<Real>::NA();
  return 0.;
#endif
}

int Binomial::rand(int n, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return (int)R::rbinom(n, prob);
#else
  return 0;
#endif
}

Real Binomial::pdf(int x, int n, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::dbinom(x, (double)n, prob, false);
#else
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_pdf_raw(x, n, prob);
#endif
}

Real Binomial::lpdf(int x, int n, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::dbinom((double)x, (double)n, prob, true);
#else
  // trivial cases
  if (Arithmetic<Real>::isNA(x)) return Arithmetic<Real>::NA();
  // compute result
  return Funct::binomial_lpdf_raw(x, n, prob);
#endif
}
Real Binomial::cdf(Real const& t, int n, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return R::pbinom(t, (double)n, prob, true, false);
#else
  return 0.;
#endif
}

int Binomial::icdf(Real const& p, int n, Real const& prob)
{
#ifdef IS_RTKPP_LIB
  return (int)R::qbinom(p, (double)n, prob, true, false);
#else
  return 0.;
#endif
}


} // namespace Law

} // namespace STK

