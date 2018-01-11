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

/*    Project: MixAll
 * created on: Jan 11, 2018
 *     Author: iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file init-MixAll.c
 *  @brief In this file we implement the R
 **/
// RegisteringDynamic Symbols

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h> // for NULL

// all MixAll method that can be used from R are defined there
#include "learnMixAll.h"
#include "clusterMixAll.h"

// declare functions
static const R_CallMethodDef callMethods[]  =
{
  {"clusterMixture", (DL_FUNC) &clusterMixture, 6},
  {"clusterMixedData", (DL_FUNC) &clusterMixedData, 5},
  {"clusterKernelMixture", (DL_FUNC) &clusterKernelMixture, 6},
  {"cluserKernelCompute", (DL_FUNC) &clusterKernelCompute, 3},
  {"learnMixture", (DL_FUNC) &learnMixture, 6},
  {"learnMixedData", (DL_FUNC) &learnMixedData, 5},
  {"learnKernelMixture", (DL_FUNC) &learnKernelMixture, 6},
  {NULL}
};


void R_init_myRoutines(DllInfo *info)
{
  /* Register the .Call routines.
  No .C  .Fortran() or .External() routines,
  so pass those arrays as NULL.
  */
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
