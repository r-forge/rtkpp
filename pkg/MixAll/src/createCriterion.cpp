/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public
    License along with this program; if not, write to the
    Free Software Foundation, Inc.,
    59 Temple Place,
    Suite 330,
    Boston, MA 02111-1307
    USA

    Contact : S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file createCriterion.cpp
 *  @brief In this file we create the given criterion
 **/


#include "../inst/projects/MixAll/MixAll_Util.h"
#include "../inst/projects/Clustering/include/STK_Clust_Util.h"

/* @return a pointer on the class computing the criterion */
STK::IMixtureCriterion* createCriterion( std::string const& criterion)
{
  STK::Clust::criterionType type = STK::Clust::stringToCriterion(criterion);
  return STK::Clust::createCriterion(type);
}


