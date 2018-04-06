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

/** @file STK_KernelMixtureManager.h
 *  @brief In this file we define the KernelMixtureManager class.
 **/


#ifndef STK_KERNELMIXTUREMANAGER_H
#define STK_KERNELMIXTUREMANAGER_H

#include "STK_KmmBridge.h"
#include "STK_KernelHandler.h"
#include "../STK_IMixtureManager.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the MixtureComposer class.
 *
 *  It allows to handle all the creation and initialization stuff needed by the
 *  (bridged) mixture models of the stkpp library.
 */
class KernelMixtureManager: public IMixtureManager<KernelHandler>
{
  public:
    typedef IMixtureManager<KernelHandler> Base;

    typedef KmmBridge< Clust::Kmm_sk_, CSquareX> KmmBridge_sk;
    typedef KmmBridge< Clust::Kmm_s_, CSquareX> KmmBridge_s;

    /** Default constructor, need an instance of a KernelHandler.  */
    KernelMixtureManager(KernelHandler const& handler);
    /** destructor */
    virtual ~KernelMixtureManager();

    /** set the dimension of the kernel mixture model */
    void setDim(IMixture* p_mixture, Real const& dim) const;
    /** set the dimension of the kernel mixture model */
    template<class Vector>
    void setDim(IMixture* p_mixture, ExprBase<Vector> const& dim) const;
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array to return with the parameters
     **/
    virtual void getParameters(IMixture* p_mixture, ArrayXX& param) const;
    /** set the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param param the array with the parameters to set
     **/
    virtual void setParameters(IMixture* p_mixture, ArrayXX const& param) const;

  protected:
    /** create a concrete mixture from its string name and initialize it.
     *  @param modelName, idData name of the model and id of the data
     *  @param nbCluster number of cluster of the model
     **/
    virtual IMixture* createMixtureImpl(String const& modelName, String const& idData, int nbCluster);

  private:
    /** create a concrete mixture and initialize it.
     *  @param idModel,idData Id of the model and the data
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixtureImpl(Clust::Mixture idModel, String const& idData, int nbCluster);
};

/* set the dimension of the kernel mixture model */
template<class Vector>
void KernelMixtureManager::setDim(IMixture* p_mixture, ExprBase<Vector> const& dim) const
{
  STK_STATIC_ASSERT_ONE_DIMENSION_ONLY(Vector);
  if (!p_mixture) return;
  Clust::Mixture idModel = getIdModel(p_mixture->idData());
  // up-cast... (Yes it's bad....;)...)
  switch (idModel)
  {
    // Kernel models
    case Clust::Kmm_sk_:
    { static_cast<KmmBridge_sk*>(p_mixture)->setDim(dim.asDerived());}
    break;
    case Clust::Kmm_s_:
    { static_cast<KmmBridge_s*>(p_mixture)->setDim(dim.asDerived());}
    break;
    default: // idModel is not implemented
    break;
  }
}


} // namespace STK

#endif /* STK_KERNELMIXTUREMANAGER_H */
