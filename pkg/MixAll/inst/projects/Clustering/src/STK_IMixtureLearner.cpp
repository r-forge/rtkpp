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
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IMixtureLearner.cpp
 *  @brief In this file we implement the interface base class for learners.
 **/


#include "../include/STK_IMixtureLearner.h"

namespace STK
{

/* Constructor.
 *  @param nbCluster,nbSample number of clusters and samples
 **/
IMixtureLearner::IMixtureLearner( int nbSample, int nbCluster)
                                : IMixtureStatModel(nbSample, nbCluster)
                                , state_(Clust::modelCreated_)
{}
/* copy constructor.
 *  @param model the model to clone
 **/
IMixtureLearner::IMixtureLearner( IMixtureLearner const& model)
                                : IMixtureStatModel(model), state_(model.state_)
{}
/* destructor */
IMixtureLearner::~IMixtureLearner() {}

} // namespace STK
