/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2018  Serge Iovleff, Université Lille 1, Inria

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

/*    Project: stkpp::MixAll
 * created on: Mar 17, 2018
 *     Author: iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file IClusterPredictor.cpp
 *  @brief In this file we implement the IClusterPredictor class
 **/


#include "../inst/projects/MixAll/IClusterPredictor.h"

namespace STK
{
IClusterPredictor::IClusterPredictor( Rcpp::S4 s4_model,  Rcpp::S4 s4_clusterPredict)
                                  : ILauncherBase(s4_model)
                                  , s4_clusterPredict_(s4_clusterPredict)
                                  , s4_algo_(s4_clusterPredict_.slot("algo"))
                                  , facade_()
                                  , p_algo_(createAlgo())
{}
IClusterPredictor::~IClusterPredictor()
{ if (p_algo_) delete p_algo_;}

/* utility function creating STK algorithm from R algorithm */
IMixtureAlgoPredict* IClusterPredictor::createAlgo()
{
   std::string algoName = s4_algo_.slot("algo");
   int nbIterBurn = s4_algo_.slot("nbIterBurn");
   int nbIterLong = s4_algo_.slot("nbIterLong");
   double epsilon = s4_algo_.slot("epsilon");
   Clust::algoPredictType algoType = Clust::stringToPredictAlgo(algoName);
   return Clust::createPredictAlgo(algoType,nbIterBurn,nbIterLong,epsilon);
}


} // namespace SK
