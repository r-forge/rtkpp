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

#include "STK_KernelGaussianBridge.h"
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

    typedef KernelGaussianBridge< Clust::KernelGaussian_sk_, CSquareX> KernelGaussianBridge_sk;
    typedef KernelGaussianBridge< Clust::KernelGaussian_s_, CSquareX> KernelGaussianBridge_s;

    /** Default constructor, need an instance of a KernelHandler.  */
    KernelMixtureManager(KernelHandler const& handler);
    /** destructor */
    virtual ~KernelMixtureManager();

    /** set the kernel of the mixture model
     * @param p_mixture pointer on the mixture
     * @param p_kernel pointer on the kernel
     * */
    void setKernel(IMixture* p_mixture, Kernel::IKernel const* p_kernel);
    /** set the dimension of the kernel mixture model */
    void setDim(IMixture* p_mixture, Real const& dim);
    /** set the dimension of the kernel mixture model */
    void setDim(IMixture* p_mixture, CPointX const& dim);
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

} // namespace STK

#endif /* STK_KERNELMIXTUREMANAGER_H */
