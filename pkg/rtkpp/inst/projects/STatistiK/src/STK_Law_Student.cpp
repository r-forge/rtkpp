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
 * Purpose:  Student probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Student.cpp
 *  @brief In this file we implement the Student probability distribution.
 **/

#include "../include/STK_Law_Student.h"

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
Student::Student( int df) : Base(_T("Student")), df_(df) {}
/* destructor */
Student::~Student() {}
/* @return a pseudo Student random variate. */
Real Student::rand() const
{
#ifdef IS_RTKPP_LIB
  return R::rt(df_);
#else
  return 0;
#endif
}
/* @return the value of the pdf
 *  @param x a positive real value
 **/
Real Student::pdf(Real const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dt(x, df_, false);
#else
  return 0;
#endif
}
/* @return the value of the log-pdf
 *  @param x a positive real value
 **/
Real Student::lpdf(Real const& x) const
{
#ifdef IS_RTKPP_LIB
  return R::dt(x, df_, true);
#else
  return 0;
#endif
}
/* @return the cumulative distribution function
 *  @param t a positive real value
 **/
Real Student::cdf(Real const& t) const
{
#ifdef IS_RTKPP_LIB
  return R::pt(t, df_, true, false);
#else
  return 0;
#endif
}
/* @return the inverse cumulative distribution function
 *  @param p a probability number
 **/
Real Student::icdf(Real const& p) const
{
#ifdef IS_RTKPP_LIB
  return R::qt(p, df_, true, false);
#else
  return 0;
#endif
}

/* @return a pseudo Student random variate with the specified parameters.
 *  @param df degree of freedom parameter
 **/
Real Student::rand( int df)
{
#ifdef IS_RTKPP_LIB
  return R::rt(df);
#else
  return 0;
#endif
}
/* @return the value of the pdf
 *  @param x a positive real value
 *  @param df degree of freedom parameter
 **/
Real Student::pdf(Real const& x, int df)
{
#ifdef IS_RTKPP_LIB
  return R::dt(x, df, false);
#else
  return 0;
#endif
}
/* @return the value of the log-pdf
 *  @param x a positive real value
 *  @param df degree of freedom parameter
 **/
Real Student::lpdf(Real const& x, int df)
{
#ifdef IS_RTKPP_LIB
  return R::dt(x, df, true);
#else
  return 0;
#endif
}
/* @return the cumulative distribution function
 *  @param t a positive real value
 *  @param df degree of freedom parameter
 **/
Real Student::cdf(Real const& t, int df)
{
#ifdef IS_RTKPP_LIB
  return R::pt(t, df, true, false);
#else
  return 0;
#endif
}
/* @return the inverse cumulative distribution function
 *  @param p a probability number
 *  @param df degree of freedom parameter
 **/
Real Student::icdf(Real const& p, int df)
{
#ifdef IS_RTKPP_LIB
  return R::qt(p, df, true, false);
#else
  return 0;
#endif
}

} // namespace Law

} // namespace STK
