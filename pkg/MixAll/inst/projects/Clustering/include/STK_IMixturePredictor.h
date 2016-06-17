/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff

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
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixturePredictor.h
 *  @brief In this file we define the abstract base class for predictors.
 **/

#ifndef STK_IMIXTUREPREDICTOR_H
#define STK_IMIXTUREPREDICTOR_H

#include "STK_Clust_Util.h"

#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <Arrays/include/STK_CArray.h>

namespace STK
{

/** @ingroup Clustering
 *  @brief Base class for Predictor of a Mixture (composed) model.
 *
 * The pure virtual function to implement in derived class are
 * @code
 *   virtual IMixturePredictor* create() const = 0;
 *   virtual IMixturePredictor* clone() const = 0;
 *   virtual Real lnComponentProbability(int i, int k) = 0;
 * @endcode
 *
 * @sa IMixtureComposer
 */
class IMixturePredictor
{
  protected:
    /** Constructor.
     **/
    IMixturePredictor( );

    /** copy constructor. If the pointer on the mixture parameters are not zero
     *  then they are cloned.
     *  @param model the model to clone
     **/
    IMixturePredictor( IMixturePredictor const& model);

  public:
    /** destructor */
    ~IMixturePredictor();

    /** @return a pointer on the proportions of each mixtures */
    inline CPointX const* p_pk() const { return p_pk_;};
    /** @return a pointer on the sum of the columns of tik = estimated number of individuals */
    inline CPointX const* p_nk() const { return p_nk_;};
    /** @return a pointer on the the tik probabilities */
    inline CArrayXX const* p_tik() const { return p_tik_;};
    /** @return a pointer on the zi class labels */
    inline CVectorXi const* p_zi() const { return p_zi_;};

  protected:
    /** The proportions of each mixtures */
    CPointX* p_pk_;
    /** The tik probabilities */
    CArrayXX* p_tik_;
    /** The sum of the columns of tik_ */
    CPointX* p_nk_;
    /** The zi class label */
    CVectorXi* p_zi_;
};

} // namespace STK

#endif /* STK_IMIXTUREPREDICTOR_H */

