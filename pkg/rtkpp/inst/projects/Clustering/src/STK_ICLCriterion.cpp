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
 * Project:  stkpp::Clustering
 * created on: 19 sept. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_ICLCriterion.cpp
 *  @brief In this file we implement the run method of the ICLCriterion class.
 **/

#include "../include/STK_ICLCriterion.h"

namespace STK
{

/* implementation of the virtual method run */
bool ICLCriterion::run()
{
  if (!p_composer_)
  { msg_error_ = STKERROR_NO_ARG(BICCriterion::run,p_model_ is not set);
    return false;
  }
  value_  = - 2.*p_composer_->lnLikelihood()
            + p_composer_->nbFreeParameter()*p_composer_->lnNbSample()
           -  2.*(p_composer_->tik()*p_composer_->tik().safe(1e-15).log()).sum();
  return true;
}

} // namespace STK


