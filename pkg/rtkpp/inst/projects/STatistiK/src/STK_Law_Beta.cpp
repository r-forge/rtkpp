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

/** @file STK_Law_Bernoulli.cpp
 *  @brief In this file we implement the Beta class.
 **/

#include "../include/STK_Law_Beta.h"
#include "../include/STK_Law_Gamma.h"
#include "Analysis/include/STK_Funct_betaRatio.h"

#ifdef IS_RTKPP_LIB
#include <Rcpp.h>
#endif
//
namespace STK
{
namespace Law
{
/* Ctor
 */
Beta::Beta( Real const& alpha, Real const& beta)
          : IUnivLaw<Real>(String(_T("Beta")))
          , alpha_(alpha)
          , beta_(beta)
{
  // check parameters
  if ( !isFinite(alpha) || !isFinite(beta) || alpha <= 0.0 || beta <= 0.0)
    STKDOMAIN_ERROR_2ARG("Beta::Beta",alpha,beta,"argument error");
}

/* Dtor
 */
Beta::~Beta() {}

/*  Generate a pseudo Beta random variate.
 */
Real Beta::rand() const
{
#ifdef IS_RTKPP_LIB
  return R::rbeta(alpha_, beta_);
#else
  Real g1 = Law::Gamma::rand(alpha_, 1.);
  return g1/(g1+Law::Gamma::rand(beta_, 1.));
#endif
}

/*
 *  Give the value of the pdf at x.
 */
Real Beta::pdf( Real const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dbeta(x,alpha_, beta_, false);
#else
  // trivial cases
  if (!Arithmetic<Real>::isFinite(x)||(x<0.)||(x>1)) return 0.0;
  // compute result
  return 0.;
#endif
}

/*
 * Give the value of the log-pdf at x.
 */
Real Beta::lpdf( Real const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dbeta(x,alpha_, beta_, true);
#else
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // check parameter
  if (Arithmetic<Real>::isInfinite(x))
    return -Arithmetic<Real>::infinity();
  // compute result
  return 0.;
#endif
}

/*
 * The cumulative distribution function at t.
 */
Real Beta::cdf( Real const& t) const
{
#ifdef IS_RTKPP_LIB
  return R::pbeta(t, alpha_, beta_, true, false);
#else
  // check NA value
  if (isNA(t)) return Arithmetic<Real>::NA();
  // compute result
  return (Arithmetic<Real>::isInfinite(t)) ? (t < 0.) ? 0.0 : 1.0
                                           :  Funct::betaRatio(alpha_,beta_,t, false);
#endif
}
    
/*
 * The inverse cumulative distribution function at p.
 */
Real Beta::icdf( Real const& p) const
{
#ifdef IS_RTKPP_LIB
  return R::qbeta(p , alpha_, beta_, true, false);
#else
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Beta::icdf,p,argument outside [0;1]);
  // trivial cases
  if (p == 0.) return 0.;
  if (p == 1.) return 1.;
  // result 
  return 0.;
#endif
}

/*  Generate a pseudo Beta random variate with the specified parameters.
 *  (static)
 */
Real Beta::rand( Real const& alpha, Real const& beta)
{
#ifdef IS_RTKPP_LIB
  return R::rbeta(alpha, beta);
#else
  Real g1 = Law::Gamma::rand(alpha, 1.);
  return g1/(g1+Law::Gamma::rand(beta, 1.));
#endif
}

Real Beta::pdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef IS_RTKPP_LIB
  return R::dbeta(x,alpha, beta, false);
#else
  // trivial cases
  if (!Arithmetic<Real>::isFinite(x)||(x<0.)||(x>1)) return 0.0;
  // compute result
  return 0.;
#endif
}

Real Beta::lpdf(const Real& x, const Real& alpha, const Real& beta)
{
#ifdef IS_RTKPP_LIB
  return R::dbeta(x,alpha, beta, true);
#else
  // check NA value
  if (isNA(x)) return Arithmetic<Real>::NA();
  // check parameter
  if (Arithmetic<Real>::isInfinite(x))
    return -Arithmetic<Real>::infinity();
  // compute result
  return 0.;
#endif
}

Real Beta::cdf(const Real& t, const Real& alpha, const Real& beta)
{
#ifdef IS_RTKPP_LIB
  return R::pbeta(t, alpha, beta, true, false);
#else
  // check NA value
  if (isNA(t)) return Arithmetic<Real>::NA();
  // compute result
  return (Arithmetic<Real>::isInfinite(t)) ? (t < 0.) ? 0.0 : 1.0
                                           :  Funct::betaRatio(alpha,beta,t, false);
#endif
}

Real Beta::icdf(const Real& p, const Real& alpha, const Real& beta)
{
#ifdef IS_RTKPP_LIB
  return R::qbeta(p , alpha, beta, true, false);
#else
  // check parameter
  if ((p > 1.) || (p < 0.))
    STKDOMAIN_ERROR_1ARG(Beta::icdf,p,argument outside [0;1]);
  // trivial cases
  if (p == 0.) return 0.;
  if (p == 1.) return 1.;
  // result
  return 0.;
#endif
}

} // namespace Law

} // namespace STK

