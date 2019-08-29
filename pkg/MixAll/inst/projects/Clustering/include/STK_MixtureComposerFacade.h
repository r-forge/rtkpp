/*--------------------------------------------------------------------*/
/*  Copyright (C) 2004-2015  Serge Iovleff, University Lille 1, Inria

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
 * created on: 11 Mars 2018
 * Author:   Iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureComposerFacade.h
 *  @brief In this file we define a facade class for p_composers.
 **/


#ifndef STK_MIXTURECOMPOSERFACADE_H
#define STK_MIXTURECOMPOSERFACADE_H

#include "./CategoricalModels/STK_CategoricalMixtureManager.h"
#include "./DiagGaussianModels/STK_DiagGaussianMixtureManager.h"
#include "./GammaModels/STK_GammaMixtureManager.h"
#include "./KernelModels/STK_KernelMixtureManager.h"
#include "./PoissonModels/STK_PoissonMixtureManager.h"

namespace STK
{

/** The MixtureComposerFacade allows to interact with the p_composer or learner
 *  for estimating a mixture model with less effort.
 *  This class stores publicly an instance of a DataHandler and a KernelHandler.
 **/
template<class DataHandler_>
class MixtureComposerFacade
{
  public:
    /** constructor.
     *  @param p_composer p_composer to use
     **/
    MixtureComposerFacade( IMixtureComposer* p_composer = 0)
                         : p_composer_(p_composer)
                         , handler_()
                         , kmmManager_(kerHandler_)
                         , kerHandler_()
                         , diagGaussianManager_(handler_)
                         , poissonManager_(handler_)
                         , gammaManager_(handler_)
                         , categoricalManager_(handler_)
    {}
    /** destructor. */
    ~MixtureComposerFacade() { if (p_composer_) delete p_composer_;}

    // not const getter. Get access to the composer, data handler and kernel manager
    /** @return pointer on the main p_composer (not const) */
    inline IMixtureComposer*& p_composer() { return p_composer_;};
    /** data handler for model based mixture models */
    inline DataHandler_& handler() { return handler_;};
    /** @return Kernel Mixture Manager (not const) */
    inline KernelMixtureManager& kmmManager() { return kmmManager_;}

    // useful
    /** add an instance of a kernel to the kernel handler */
    bool addKernel(Kernel::IKernel* p_kernel,  String const& idData, String const& idModel);
    /** create composer */
    void createComposer(int nbSample, int nbCluster);
    /** create fixed proportion composer */
    void createFixedPropComposer(int nbSample, int nbCluster);
    /** delete composer */
    void deleteComposer();
    /** create all mixtures using managers */
    void createMixtures();

    // const getters
    /** @return Kernel Handler mixtures manager*/
    inline KernelHandler const& kerHandler() const { return kerHandler_;}
    /** @return DiagGaussian mixtures manager*/
    inline DiagGaussianMixtureManager<DataHandler_> const& diagGaussianManager() const
    { return diagGaussianManager_;}
    /** @return Poisson mixture models manager */
    inline PoissonMixtureManager<DataHandler_> const& poissonManager() const
    { return poissonManager_;}
    /** gamma mixture models manager */
    inline GammaMixtureManager<DataHandler_> const& gammaManager() const
    { return gammaManager_;}
    /** categorical mixture models manager */
    inline CategoricalMixtureManager<DataHandler_> const& categoricalManager() const
    { return categoricalManager_;}

    /** get parameters */
    bool getParameters( String const& idData, ArrayXX& param) const;
    /** get diagonal Gaussian parameters */
    void getDiagGaussianParameters( String const& idData, ArrayXX& param) const;
    /** get Poisson parameters */
    void getPoissonParameters( String const& idData, ArrayXX& param) const;
    /** get gamma parameters */
    void getGammaParameters( String const& idData, ArrayXX& param) const;
    /** get categorical parameters */
    void getCategoricalParameters( String const& idData, ArrayXX& param) const;
    /** get Kernel Mixture Model parameters */
    void getKmmParameters( String const& idData, ArrayXX& param) const;

    // setter
    /** set composer */
    void setComposer(IMixtureComposer* p_composer);
    /** set parameters of a component */
    bool setParameters( String const& idData, ArrayXX const& param);
    /** set diagonal Gaussian parameters */
    void setDiagGaussianParameters( String const& idData, ArrayXX const& param);
    /** set Poisson parameters */
    void setPoissonParameters( String const& idData, ArrayXX const& param);
    /** set gamma parameters */
    void setGammaParameters( String const& idData, ArrayXX const& param);
    /** set categorical parameters */
    void setCategoricalParameters( String const& idData, ArrayXX const& param);
    /** set Kernel Mixture Model parameters */
    void setKmmParameters( String const& idData, ArrayXX const& param);

    // template setters
    /** set the mixture parameters using an array of posterior probabilities.
     *  Proportions, numbers in each class and class labels are computed
     *  using these posterior probabilities.
     *  @param tik posterior class probabilities
     **/
    template<class Array>
    void setMixtureParameters( Array const& tik);
    /** set the mixture parameters giving the posterior probabilities and
     *  the proportions.
     *  Numbers in each class and class labels are computed using the
     *  posterior probabilities.
     *  @param tik posterior class probabilities
     *  @param pk prior class proportion
     **/
    template<class Array, class RowVector>
    void setMixtureParameters( Array const& tik, RowVector const& pk);
    /** Set proportions of each classes
     *  @param pk prior class proportion
     **/
    template<class RowVector>
    void setProportions( RowVector const& pk);

  protected:
    /** pointer on the main p_composer */
    IMixtureComposer* p_composer_;
    /** data handler for model based mixture models */
    DataHandler_ handler_;
    /** Kernel Mixture Manager */
    KernelMixtureManager kmmManager_;
    /** handler for kernel mixture models */
    KernelHandler kerHandler_;
    /** diagonal Gaussian mixture models manager */
    DiagGaussianMixtureManager<DataHandler_> diagGaussianManager_;
    /** Poisson mixture models manager */
    PoissonMixtureManager<DataHandler_> poissonManager_;
    /** gamma mixture models manager */
    GammaMixtureManager<DataHandler_> gammaManager_;
    /** categorical mixture models manager */
    CategoricalMixtureManager<DataHandler_> categoricalManager_;
};

/* add an instance of a kernel to the kernel handler */
template<class DataHandler_>
bool MixtureComposerFacade<DataHandler_>::addKernel(Kernel::IKernel* p_kernel,  String const& idData, String const& idModel)
{ return kerHandler_.addKernel(p_kernel, idData, idModel);}
/* create composer */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::createComposer(int nbSample, int nbCluster)
{ p_composer_ = new MixtureComposer(nbSample, nbCluster);}
/* create composer */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::createFixedPropComposer(int nbSample, int nbCluster)
{ p_composer_ = new MixtureComposerFixedProp(nbSample, nbCluster);}
/* create composer */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::deleteComposer()
{
  if (p_composer_) { delete p_composer_ ; p_composer_ = 0;}
}
/* create the mixtures in the given learner */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::createMixtures()
{
  p_composer_->createMixture(diagGaussianManager_);
  p_composer_->createMixture(poissonManager_);
  p_composer_->createMixture(gammaManager_);
  p_composer_->createMixture(categoricalManager_);
  p_composer_->createMixture(kmmManager_);
}

/* fill param with the parameters */
template<class DataHandler_>
bool MixtureComposerFacade<DataHandler_>::getParameters( String const& idData, ArrayXX& param) const
{
  // get idModel from idData
  String idModel;
  if(!handler_.getIdModelName(idData, idModel)) { return false;} // should not happened

  // get mixture id
  Clust::Mixture mix = Clust::stringToMixture(idModel);
  if (mix == Clust::unknown_mixture_)
  {
    bool prop;
    mix = Clust::stringToMixture(idModel, prop);
  }

  // get mixture class (gamma, Gaussian, Poisson, etc.) and parameters for this class of mixture
  switch (Clust::mixtureToMixtureClass(mix))
  {
    case Clust::DiagGaussian_:
      getDiagGaussianParameters( idData, param);
      return true;
      break;
    case Clust::Poisson_:
      getPoissonParameters( idData, param);
      return true;
      break;
    case Clust::Gamma_:
      getGammaParameters( idData, param);
      return true;
      break;
    case Clust::Categorical_:
      getCategoricalParameters( idData, param);
      return true;
      break;
    case Clust::Kmm_:
      getKmmParameters( idData, param);
      return true;
      break;
    case Clust::unknown_mixture_class_:
      return false;
      break;
    default:
      return false;
      break;
  }
  return false; // avoid compiler warning
}

/* set composer */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::setComposer(IMixtureComposer* p_composer)
{ p_composer_ = p_composer;}

/* set model parameters with param */
template<class DataHandler_>
bool MixtureComposerFacade<DataHandler_>::setParameters( String const& idData, ArrayXX const& param)
{
    // get idModel from idData
    String idModel;
    if(!handler_.getIdModelName(idData, idModel)) { return false;}
    // get mixture id
    Clust::Mixture mix = Clust::stringToMixture(idModel);
    if (mix == Clust::unknown_mixture_) { return false;}

    // get mixture class (gamma, Gaussian, Poisson, etc.) and parameters for this class of mixture
    switch (Clust::mixtureToMixtureClass(mix))
    {
      case Clust::DiagGaussian_:
        setDiagGaussianParameters( idData, param);
        return true;
        break;
      case Clust::Poisson_:
        setPoissonParameters( idData, param);
        return true;
        break;
      case Clust::Gamma_:
        setGammaParameters( idData, param);
        return true;
        break;
      case Clust::Categorical_:
        setCategoricalParameters( idData, param);
        return true;
        break;
      case Clust::Kmm_:
        setKmmParameters( idData, param);
        return true;
        break;
      case Clust::unknown_mixture_class_:
        return false;
        break;
      default:
        return false;
        break;
    }
    return true; // avoid compiler warning
}
/* get the diagonal Gaussian parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::getDiagGaussianParameters( String const& idData, ArrayXX& param) const
{ p_composer_->getParameters(diagGaussianManager_,idData, param);}
/* get the Poisson parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::getPoissonParameters( String const& idData, ArrayXX& param) const
{ p_composer_->getParameters(poissonManager_,idData, param);}
/* get the gamma parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::getGammaParameters( String const& idData, ArrayXX& param) const
{ p_composer_->getParameters(gammaManager_,idData, param);}
/* get the Categorical parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::getCategoricalParameters( String const& idData, ArrayXX& param) const
{ p_composer_->getParameters(categoricalManager_,idData, param);}
/* get the kernel parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::getKmmParameters( String const& idData, ArrayXX& param) const
{ p_composer_->getParameters(kmmManager_,idData, param);}


// setters
/* set the diagonal Gaussian parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::setDiagGaussianParameters( String const& idData, ArrayXX const& param)
{ p_composer_->setParameters(diagGaussianManager_,idData, param);}
/* set the diagonal Gaussian parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::setPoissonParameters( String const& idData, ArrayXX const& param)
{ p_composer_->setParameters(poissonManager_,idData, param);}
/* set the gamma parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::setGammaParameters( String const& idData, ArrayXX const& param)
{ p_composer_->setParameters(gammaManager_,idData, param);}
/* set the Categorical parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::setCategoricalParameters( String const& idData, ArrayXX const& param)
{ p_composer_->setParameters(categoricalManager_,idData, param);}
/* set the Categorical parameters */
template<class DataHandler_>
void MixtureComposerFacade<DataHandler_>::setKmmParameters( String const& idData, ArrayXX const& param)
{ p_composer_->setParameters(kmmManager_,idData, param);}

template<class DataHandler_>
template<class Array>
void MixtureComposerFacade<DataHandler_>::setMixtureParameters(Array const& tik)
{ p_composer_->setMixtureParameters(tik);}

template<class DataHandler_>
template<class Array, class RowVector>
void MixtureComposerFacade<DataHandler_>::setMixtureParameters( Array const& tik, RowVector const& pk)
{ p_composer_->setMixtureParameters(tik, pk);}

template<class DataHandler_>
template<class RowVector>
void MixtureComposerFacade<DataHandler_>::setProportions( RowVector const& pk)
{ p_composer_->setProportions(pk);}


} // namespace STK

#endif /* STK_MIXTURECOMPOSERFACADE_H */
