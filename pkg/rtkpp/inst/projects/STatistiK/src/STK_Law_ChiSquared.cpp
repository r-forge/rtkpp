/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2008  Serge Iovleff

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
 * Purpose:  ChiSquared probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_ChiSquared.cpp
 *  @brief In this file we implement the ChiSquared probability distribution.
 **/

#include "../include/STK_Law_ChiSquared.h"

#ifdef IS_RTKPP_LIB
#include <Rcpp.h>
#endif

namespace STK
{
namespace Law
{
/* Default constructor.
 *  @param df degree of freedom parameter
 **/
ChiSquared::ChiSquared(int df) : Base(_T("Chi-squared")), df_(df)
{
  if (df<=0) STKDOMAIN_ERROR_1ARG(ChiSquared::ChiSquared,df,df must be > 0);
}
/* destructor */
ChiSquared::~ChiSquared() {}
/* @return a pseudo ChiSquared random variate. */
Real ChiSquared::rand() const
{
#ifdef IS_RTKPP_LIB
  return R::rchisq(df_);
#else
  return 0.;
#endif
}
/* @return the value of the pdf
 *  @param x a positive real value
 **/
Real ChiSquared::pdf(const Real& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dchisq(x, df_, false);
#else
  return 0.;
#endif
}

Real ChiSquared::lpdf(const Real& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dchisq(x, df_, true);
#else
  return 0.;
#endif
}

/* @return the cumulative distribution function
 *  @param t a positive real value
 **/
Real ChiSquared::cdf(const Real& t) const
{
#ifdef IS_RTKPP_LIB
  return R::pchisq(t, df_, true, false);
#else
  return 0.;
#endif
}

/* @return the inverse cumulative distribution function
 *  @param p a probability number
 **/
Real ChiSquared::icdf(const Real& p) const
{
#ifdef IS_RTKPP_LIB
  return R::qchisq(p, df_, true, false);
#else
  return 0.;
#endif
}

/* @return a pseudo ChiSquared random variate with the specified parameters.
 *  @param df degree of freedom parameter
 **/
Real ChiSquared::rand(int df)
{
#ifdef IS_RTKPP_LIB
  return R::rchisq(df);
#else
  return 0.;
#endif
}

/* @return the value of the pdf
 *  @param x a positive real value
 *  @param df degree of freedom parameter
 **/
Real ChiSquared::pdf(const Real& x, int df)
{
#ifdef IS_RTKPP_LIB
  return R::dchisq(x, df, false);
#else
  return 0.;
#endif
}

/* @return the value of the log-pdf
 *  @param x a positive real value
 *  @param df degree of freedom parameter
 **/
Real ChiSquared::lpdf(const Real& x, int df)
{
#ifdef IS_RTKPP_LIB
  return R::dchisq(x, df, true);
#else
  return 0.;
#endif
}

/* @return the cumulative distribution function
 *  @param t a positive real value
 *  @param df degree of freedom parameter
 **/
Real ChiSquared::cdf(const Real& t, int df)
{
#ifdef IS_RTKPP_LIB
  return R::pchisq(t, df, true, false);
#else
  return 0.;
#endif
}

/* @return the inverse cumulative distribution function
 *  @param p a probability number
 *  @param df degree of freedom parameter
 **/
Real ChiSquared::icdf(const Real& p, int df)
{
#ifdef IS_RTKPP_LIB
  return R::qchisq(p, df, true, false);
#else
  return 0.;
#endif
}

} // namespace Law

} // namespace STK

