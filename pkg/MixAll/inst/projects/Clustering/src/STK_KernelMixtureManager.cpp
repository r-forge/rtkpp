/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016  Serge Iovleff, Université Lille 1, Inria

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
 * created on: 15 mars 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_KernelMixtureManager.cpp
 *  @brief In this file we implement the KernelMixtureManager class.
 **/


#include "../include/KernelModels/STK_KernelMixtureManager.h"

namespace STK
{

/* Default constructor, need an instance of a DataHandler.  */
KernelMixtureManager::KernelMixtureManager(KernelHandler const& handler): Base(&handler) {}
/* destructor */
KernelMixtureManager::~KernelMixtureManager() {}

/* set the kernel of the mixture model
 * @param p_mixture pointer on the mixture
 * @param p_kernel pointer on the kernel
 * */
void KernelMixtureManager::setKernel(IMixture* p_mixture, Kernel::IKernel const* p_kernel)
{
  if (!p_mixture) return;
  Clust::Mixture idModel = getIdModel(p_mixture->idData());
  // up-cast... (Yes it's bad....;)...)
  switch (idModel)
  {
    // Kernel models
    case Clust::KernelGaussian_sk_:
    { static_cast<KernelGaussianBridge_sk*>(p_mixture)->setKernel(p_kernel);}
    break;
    case Clust::KernelGaussian_s_:
    { static_cast<KernelGaussianBridge_s*>(p_mixture)->setKernel(p_kernel);}
    break;
    default: // idModel is not implemented
    break;
  }
}
/* set the dimension of the kernel mixture model */
void KernelMixtureManager::setDim(IMixture* p_mixture, Real const& dim) const
{
  if (!p_mixture) return;
  Clust::Mixture idModel = getIdModel(p_mixture->idData());
  // up-cast... (Yes it's bad....;)...)
  switch (idModel)
  {
    // Kernel models
    case Clust::KernelGaussian_sk_:
    { static_cast<KernelGaussianBridge_sk*>(p_mixture)->setDim(dim);}
    break;
    case Clust::KernelGaussian_s_:
    { static_cast<KernelGaussianBridge_s*>(p_mixture)->setDim(dim);}
    break;
    default: // idModel is not implemented
    break;
  }
}

/* set the dimension of the kernel mixture model */
void KernelMixtureManager::setDim(IMixture* p_mixture, CPointX const& dim) const
{
  if (!p_mixture) return;
  Clust::Mixture idModel = getIdModel(p_mixture->idData());
  // up-cast... (Yes it's bad....;)...)
  switch (idModel)
  {
    // Kernel models
    case Clust::KernelGaussian_sk_:
    { static_cast<KernelGaussianBridge_sk*>(p_mixture)->setDim(dim);}
    break;
    case Clust::KernelGaussian_s_:
    { static_cast<KernelGaussianBridge_s*>(p_mixture)->setDim(dim);}
    break;
    default: // idModel is not implemented
    break;
  }
}

/* get the parameters from an IMixture.
 *  @param p_mixture pointer on the mixture
 *  @param param the array to return with the parameters
 **/
void KernelMixtureManager::getParameters(IMixture* p_mixture, ArrayXX& param) const
{
  if (!p_mixture) return;
  Clust::Mixture idModel = getIdModel(p_mixture->idData());
  // up-cast... (Yes it's bad....;)...)
  switch (idModel)
  {
    // Kernel models
    case Clust::KernelGaussian_sk_:
    { static_cast<KernelGaussianBridge_sk*>(p_mixture)->getParameters(param);}
    break;
    case Clust::KernelGaussian_s_:
    { static_cast<KernelGaussianBridge_s*>(p_mixture)->getParameters(param);}
    break;
    default: // idModel is not implemented
    break;
  }
}
/* set the parameters from an IMixture.
 *  @param p_mixture pointer on the mixture
 *  @param param the array with the parameters to set
 **/
void KernelMixtureManager::setParameters(IMixture* p_mixture, ArrayXX const& param) const
{
  if (!p_mixture) return;
  Clust::Mixture idModel = getIdModel(p_mixture->idData());
  // up-cast... (Yes it's bad....;)...)
  switch (idModel)
  {
    // Kernel models
    case Clust::KernelGaussian_sk_:
    { static_cast<KernelGaussianBridge_sk*>(p_mixture)->setParameters(param);}
    break;
    case Clust::KernelGaussian_s_:
    { static_cast<KernelGaussianBridge_s*>(p_mixture)->setParameters(param);}
    break;
    default: // idModel is not implemented
    break;
  }
}

IMixture* KernelMixtureManager::createMixtureImpl(String const& modelName, String const& idData, int nbCluster)
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("KernelMixtureManager::Entering createMixtureImpl(") <<  modelName
           << _T(",") << idData << _T(",") << nbCluster << _T(")\n");
#endif
  Clust::Mixture idModel = Clust::stringToMixture(modelName);
  return createMixtureImpl(idModel, idData, nbCluster);
}
/* create a concrete mixture and initialize it.
 *  @param idModel,idData Id of the model and the data
 *  @param nbCluster number of cluster of the model
 **/
IMixture* KernelMixtureManager::createMixtureImpl(Clust::Mixture idModel, String const& idData, int nbCluster)
{
  switch (idModel)
  {
    case Clust::KernelGaussian_sk_:
    { return new KernelGaussianBridge_sk( 0, idData, nbCluster);}
    break;
    case Clust::KernelGaussian_s_:
    { return new KernelGaussianBridge_s( 0, idData, nbCluster);}
    break;
    default:
      return 0; // 0 if idModel does not exists
      break;
  }
  return 0; // 0 if idModel is not a STK++ model
}

} // namespace STK

