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
 * Purpose:  Beta probability distribution.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Law_Beta.h
 *  @brief In this file we define the Beta probability distribution.
 **/

#ifndef STK_LAW_BETA_H
#define STK_LAW_BETA_H

#include "STK_Law_IUnivLaw.h"

namespace STK
{

namespace Law
{
/** @ingroup Laws
 *  @brief Beta distribution law.
 * 
 *  In probability theory and statistics, the <em> beta distribution </em> is a
 *  family of continuous probability distributions defined on the interval [0, 1]
 *  parameterized by two positive shape parameters that appear as exponents of
 *  the random variable and control the shape of the distribution.
 *
 *  The Beta distribution, is a continuous probability distribution with pdf
 *  \f[
 *  f(x) = \frac{x^{\alpha-1} (1-x)^{\beta-1} }
 *              {{\boldsymbol \beta}(\alpha,\beta)},
 *  \quad x \in [0,1],\ \mbox{ et }\ \alpha,\beta>0
 *  \f]
 *  where \f${\boldsymbol \beta}(a,b) = \frac{\Gamma(\alpha)\Gamma(\beta)}{\Gamma(\alpha+\beta)} \f$
 *  represents the beta function.
 *
 *  @sa STK::Funct::betaRatio
**/
class Beta : public IUnivLaw<Real>
{
  public:
    /** default constructor. */
    Beta( Real const& alpha = .5, Real const& beta = .5);
    /** Dtor. */
    virtual ~Beta();
    /** @return the alpha value */
    inline Real alpha() const {return alpha_;}
    /** @return the beta value */
    inline Real beta() const { return beta_;}
    /** set the alpha value */
    inline void setAlpha(Real alpha) { alpha_ = alpha;}
    /** set the beta value */
    inline void setBeta(Real beta) {beta_ =beta;}

    /** @return a pseudo Beta random variate.
     *  This function use the Gamma::rand() random generator.
     *  TODO : implement the order statistics when a and b
     *  are small integers (Devroye p. 431). 
     *  TODO : implement the rejection method from the normal
     *  pdf when a=b (Devroye p. 434).
     *  TODO : Implement Cheng's algorithm BA for beta pdf 
     *  (Devroye p. 438)
     **/
    virtual Real rand() const;
    /** @return the value of the pdf at x. */
    virtual Real pdf( Real const& x) const;
    /** @return the value of the log-pdf at x. */
    virtual Real lpdf( Real const& x) const;
    /** @return The cumulative distribution function */
    virtual Real cdf( Real const& t) const;
    /** @return The inverse cumulative distribution */
    virtual Real icdf( Real const& p) const;

    /** @return a pseudo Beta random variate with the specified parameters. */
    static Real rand( Real const& alpha, Real const& beta);
    /** @return the value of the pdf at x. */
    static Real pdf( Real const& x, Real const& alpha, Real const& beta);
    /** @return the value of the log-pdf at x. */
    static Real lpdf( Real const& x, Real const& alpha, Real const& beta);
    /** @return The cumulative distribution function */
    static Real cdf( Real const& t, Real const& alpha, Real const& beta);
    /** @return The inverse cumulative distribution */
    static Real icdf( Real const& p, Real const& alpha, Real const& beta);

  protected:
    /** First parameter. */
    Real alpha_;
    /** Second parameter. */
    Real beta_;
};

} // namespace Law

} // namespace STK

#endif /*STK_LAWBETA_H*/
