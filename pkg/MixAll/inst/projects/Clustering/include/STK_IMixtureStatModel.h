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

/** @file STK_IMixtureStatModel.h
 *  @brief In this file we define the abstract base class for mixture statistical models.
 **/

#ifndef STK_IMIXTURESTATMODEL_H
#define STK_IMIXTURESTATMODEL_H

#include <StatModels/include/STK_IStatModelBase.h>

#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <Arrays/include/STK_CArray.h>

#include "STK_Clust_Util.h"
#include "STK_IMixture.h"
#include "STK_IMixtureManager.h"


namespace STK
{
// forward declaration
class IMixture;
/** @ingroup StatModels
 *  @brief Base class for Mixture (composed) model.
 *
 * In statistics, a mixture model is a probabilistic model for representing
 * the presence of sub-populations within an overall population, without
 * requiring that an observed data-set should identify the sub-population to
 * which an individual observation belongs. Formally a mixture model
 * corresponds to the mixture distribution that represents the probability
 * distribution of observations in the overall population. However, while
 * problems associated with "mixture distributions" relate to deriving the
 * properties of the overall population from those of the sub-populations,
 * "mixture models" are used to make statistical inferences about the
 * properties of the sub-populations given only observations on the pooled
 * population, without sub-population-identity information.
 *
 * Some ways of implementing mixture models involve steps that attribute
 * postulated sub-population-identities to individual observations (or weights
 * towards such sub-populations), in which case these can be regarded as types
 * unsupervised learning or clustering procedures. However not all inference
 * procedures involve such steps.
 *
 * The pure virtual function to implement in derived class are
 * @code
 *   virtual IMixtureComposer* create() const = 0;
 *   virtual IMixtureComposer* clone() const = 0;
 *   virtual void registerMixture(IMixture* p_mixture) = 0;
 *   virtual void releaseMixture(String const& idData) = 0;
 *   virtual bool randomInit() =0;
 *   virtual Real lnComponentProbability(int i, int k) const = 0;
 *   virtual int computeNbFreeParameters() const = 0;
 * @endcode
 *
 * The virtual function that can be re-implemented in derived class for a
 * specific behavior are:
 * @code
 *   virtual void imputationStep();
 *   virtual void samplingStep();
 *   virtual void setParametersStep();
 *   virtual void storeIntermediateResults(int iteration);
 *   virtual void releaseIntermediateResults()
 *   virtual void writeParameters(std::ostream& os) const;
 * @endcode
 *
 * @sa IMixtureComposer, IMixtureLearner
 *
 * @note the virtual method @c IMixtureStatModel::initializeStep is called in all
 * the initialization method.  Don't forget to called it in the randomInit
 * implementation.
 */
class IMixtureStatModel: public IStatModelBase
{
  protected:
    /** Constructor.
     * @param nbCluster,nbSample number of clusters and samples
     **/
    IMixtureStatModel( int nbSample, int nbCluster);
    /** copy constructor.
     *  @param model the model to clone
     **/
    IMixtureStatModel( IMixtureStatModel const& model);

  public:
    typedef std::vector<IMixture*>::const_iterator ConstMixtIterator;
    typedef std::vector<IMixture*>::iterator MixtIterator;
    /** destructor */
    virtual ~IMixtureStatModel();

    /** @return the number of cluster */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the proportions of each mixtures */
    inline CPointX const& pk() const { return pk_;};
    /** @return the tik probabilities */
    inline CArrayXX const& tik() const { return tik_;};
    /** @return the proportions of individuals */
    inline CPointX const& nk() const { return nk_;};
    /** @return the zi class label */
    inline CVectorXi const& zi() const { return zi_;};
    /** @return a constant reference on the vector of mixture */
    inline std::vector<IMixture*> const& v_mixtures() const { return v_mixtures_;}

    /** @return the computed log-likelihood of the i-th sample.
     *  @param i index of the sample
     **/
    Real computeLnLikelihood(int i) const;
    /** @return the computed likelihood of the i-th sample.
     *  @param i index of the sample
     **/
    Real computeLikelihood(int i) const;
    /** @return the computed log-likelihood. */
    Real computeLnLikelihood() const;
    /** @return the computed ICL criteria. */
    Real computeICL() const;

    /** Utility lookup function allowing to find a Mixture from its idData
     *  @param idData the id name of the data
     *  @return a pointer on the mixture, NULL if the mixture is not found
     **/
     IMixture* getMixture( String const& idData) const;
     /** @brief register a mixture to the composer.
      *  When a mixture is registered, the composer:
      *  - assign composer pointer (itself) to the mixture
      *  - add it to v_mixtures_
      *  - update the number of variables
      * @param p_mixture a pointer on the mixture
      **/
     void registerMixture(IMixture* p_mixture);
     /** @brief release a mixture from the composer.
      *  When a mixture is released, the composer remove it from v_mixtures_.
      *  @note the data set associated with the mixture is still in the manager
      *  list of data set if you use a manager in order to register it.
      *  @param idData the Id of the mixture to release.
      **/
     void releaseMixture(String const& idData);
     /** compute the number of free parameters of the model.
      *  This method is used in IMixtureComposer::initializeStep
      *  in order to give a value to IStatModelBase::nbFreeParameter_.
      *  lookup on the mixtures and sum the nbFreeParameter.
      *  @return the number of free parameters
      **/
     int computeNbFreeParameters() const;
     /** @brief compute the number of variables of the model.
      *  lookup on the mixtures and sum the nbFreeParameter.
      **/
     int computeNbVariables() const;

    // pure virtual
    /** create pattern */
    virtual IMixtureStatModel* create() const = 0;
    /** clone pattern */
    virtual IMixtureStatModel* clone() const = 0;
    /** initialize randomly the parameters of the components of the model */
    virtual void randomInit() = 0;
    /** @return the value of the probability of the i-th sample in the k-th class.
     *  @param i,k indexes of the sample and of the class
     **/
    virtual Real lnComponentProbability(int i, int k) const = 0;

    // virtual with default implementation
    /** @brief Initialize the model before at its first use.
     *  This function can be overloaded in derived class for initialization of
     *  the specific model parameters. It should be called prior to any used of
     *  the class.
     *  @sa IMixture,MixtureBridge,MixtureLearner
     **/
    virtual void initializeStep();

    /** @brief Impute the missing values.
     *  Default behavior is "do nothing".
     **/
    inline virtual void imputationStep() {}
    /** @brief Simulation of all the latent variables and/or missing data
     *  excluding class labels. Default behavior is "do nothing".
     */
    inline virtual void samplingStep() {}
    /** @brief Utility method allowing to signal to a mixture to set its parameters.
     *  It will be called once enough intermediate results have been stored.
     **/
    inline virtual void setParametersStep() {}
    /**@brief This step can be used to signal to the mixtures that they must
     * store results. This is usually called after a burn-in phase. The composer
     * store the current value of the log-Likelihood.
     **/
    inline virtual void storeIntermediateResults(int iteration) {}
    /**@brief This step can be used to signal to the mixtures that they must
     * release the stored results. This is usually called if the estimation
     * process failed.
     **/
    inline virtual void releaseIntermediateResults() {}
    /** @brief Finalize the estimation of the model.
     *  The default behavior is compute current lnLikelihood.
     **/
    inline virtual void finalizeStep()
    { setLnLikelihood(computeLnLikelihood());}
    /** write the parameters of the model in the stream os. */
    virtual void writeParameters(ostream& os) const {};

    /** Utility method allowing to create all the mixtures registered in the
     *  data handler of a mixture manager and to register them.
     *  @param manager the manager with the responsibility of the creation.
     **/
    template<class DataHandler>
    void createMixture(IMixtureManager<DataHandler>& manager);
    /** Utility method allowing to create a mixture with a given data set
     *  and register it. The Mixture Manager will find the associated model
     *  to use with this data set.
     *  @param manager the manager with the responsibility of the creation.
     *  @param idData the id name of the data to modelize.
     **/
    template<class DataHandler>
    IMixture* createMixture(IMixtureManager<DataHandler>& manager, String const& idData);
    /** Utility method allowing to release completely a mixture with its data set.
     *  The MixtureManager will find and release the associated data set.
     *  @param manager the manager with the responsibility of the release.
     *  @param idData the id name of the data to modelize.
     **/
    template<class DataHandler>
    void removeMixture(IMixtureManager<DataHandler>& manager, String const& idData);
    /** Utility method allowing to get the parameters of a specific mixture.
     *  @param manager the manager with the responsibility of the parameters
     *  @param idData the Id of the data we want the parameters
     *  @param param the structure which will receive the parameters
     **/
    template<class DataHandler, class Parameters>
    void getParameters(IMixtureManager<DataHandler> const& manager, String const& idData, Parameters& param) const
    { manager.getParameters(getMixture(idData), param);}
    /** Utility method allowing to set the parameters to a specific mixture.
     *  @param manager the manager with the responsibility of the parameters
     *  @param idData the Id of the data we want to set the parameters
     *  @param param the structure which contains the parameters
     **/
    template<class DataHandler, class Parameters>
    void setParameters(IMixtureManager<DataHandler> const& manager, String const& idData, Parameters const& param)
    { manager.setParameters(getMixture(idData), param);}

  protected:
    /** set the number of cluster of the model
     *  @param nbCluster number of cluster of the model
     * */
    inline void setNbCluster( int nbCluster) { nbCluster_ = nbCluster;}
    /** number of cluster. */
    int nbCluster_;
    /** The proportions of each mixtures */
    CPointX pk_;
    /** The tik probabilities */
    CArrayXX tik_;
    /** The sum of the columns of tik_ */
    CPointX nk_;
    /** The zi class label */
    CVectorXi zi_;

    /** vector of pointers to the mixtures components */
    std::vector<IMixture*> v_mixtures_;
};


/* Utility method allowing to create all the mixtures handled by a mixture
 * manager.
 *  @param manager the manager with the responsibility of the creation.
 **/
template<class DataHandler>
void IMixtureStatModel::createMixture(IMixtureManager<DataHandler>& manager)
{
  typedef typename DataHandlerBase<DataHandler>::InfoMap InfoMap;
  for ( typename InfoMap::const_iterator it=manager.p_handler()->info().begin(); it!=manager.p_handler()->info().end(); ++it)
  {
    IMixture* p_mixture = manager.createMixture(it->first, nbCluster());
#ifdef STK_MIXTURE_DEBUG
  if (!p_mixture)
  { stk_cout << _T("In IMixtureStatModel::createMixture(manager) failed.\n");}
#endif
    if (p_mixture) registerMixture(p_mixture);
  }
}

/* Utility method allowing to create a mixture with a given data set
 *  and register it. The Mixture Manager will find the associated model
 *  to use with this data set.
 *  @param manager the manager with the responsibility of the creation.
 *  @param idData the id name of the data to modelize.
 **/
template<class DataHandler>
IMixture* IMixtureStatModel::createMixture(IMixtureManager<DataHandler>& manager, String const& idData)
{
  IMixture* p_mixture = manager.createMixture( idData, nbCluster());
#ifdef STK_MIXTURE_DEBUG
  if (!p_mixture)
  { stk_cout << _T("In IMixtureStatModel::createMixture(manager,")<< idData << _T(") failed.\n");}
#endif
  if (p_mixture) registerMixture(p_mixture);
  return p_mixture;
}

/* Utility method allowing to release completely a mixture with its data set.
 *  The MixtureManager will find and release the associated data set.
 *  @param manager the manager with the responsibility of the release.
 *  @param idData the id name of the data to modelize.
 **/
template<class DataHandler>
void IMixtureStatModel::removeMixture(IMixtureManager<DataHandler>& manager, String const& idData)
{
  IMixture* p_mixture = getMixture(idData);
#ifdef STK_MIXTURE_DEBUG
  if (!p_mixture)
  { stk_cout << _T("In IMixtureStatModel::removeMixture(manager,")<< idData << _T(") failed.\n");}
#endif
  if (p_mixture)
  {
    releaseMixture(idData);
    manager.releaseDataBridge( idData);
  }
}

} // namespace STK

#endif /* STK_IMIXTURESTATMODEL_H */


