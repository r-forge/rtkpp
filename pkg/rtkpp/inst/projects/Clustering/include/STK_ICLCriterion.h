/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff

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
 * Project:  stkpp::
 * created on: 19 sept. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ICLCriterion.h
 *  @brief In this file we define the ICLCriterion class.
 **/


#ifndef STK_ICLCRITERION_H
#define STK_ICLCRITERION_H

#include "../../StatModels/include/STK_PenCriterion.h"
#include "STK_IMixtureComposer.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief Derived class of ICriterion for computing the ICL Criterion
 *  The ICL criteria is a penalization of the likelihood given by the formula
 *  \f[
 *  \ln{p(x|K)} \approx \mathrm{ICL} = { -2 \ln{L} + D_K \log(n) - 2 \sum_i\sum_k t_{ik}\log(t_{ik} }
 *  \f]
 *  where \f$ L \f$ represents the likelihood of the observations and \f$ D_K \f$
 *  the number of free parameters of the model.
 **/
class ICLCriterion : protected ICriterion
{
  public:
    /** Default Constructor. */
    inline ICLCriterion() : ICriterion(), p_composer_(0) {}
    /** Constructor.
     *  @param p_composer a pointer on the current composer
     **/
    inline ICLCriterion( IMixtureComposer* const p_composer)
                       : ICriterion(p_composer), p_composer_(p_composer) {}
    /** Constructor.
     *  @param composer the current composer
     **/
    inline ICLCriterion( IMixtureComposer const& composer)
                       : ICriterion(composer), p_composer_(&composer) {}
    /** copy Constructor.
     *  @param criterion the criterion to copy
     **/
    inline ICLCriterion( ICLCriterion const& criterion)
                       : ICriterion(criterion), p_composer_(criterion.p_composer_) {}
    /** virtual destructor. */
    inline virtual ~ICLCriterion() {}
    /** clone pattern */
    inline ICLCriterion* clone() const { return new ICLCriterion(*this);}
    /** @return The value of the criterion */
    inline Real const& value() const { return ICriterion::value();}
    /** @param p_composer a pointer on the current model to set */
    inline void setModel( IMixtureComposer* const p_composer)
    { ICriterion::setModel(p_composer); p_composer_ = p_composer;}
    /** @param model the current model to set */
    inline void setModel( IMixtureComposer const& model)
    { ICriterion::setModel(model); p_composer_ = &model;}
    /** implementation of the virtual method run */
    virtual bool run();

  private:
    IMixtureComposer const* p_composer_;
};


} // namespace STK

#endif /* STK_ICLCRITERION_H */
