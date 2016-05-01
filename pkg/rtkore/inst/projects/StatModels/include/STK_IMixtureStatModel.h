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

#include "STK_IStatModelBase.h"

#include <Arrays/include/STK_CArrayPoint.h>
#include <Arrays/include/STK_CArrayVector.h>
#include <Arrays/include/STK_CArray.h>

namespace STK
{

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
 *   virtual void setParameters();
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
    /** compute the number of free parameters of the model.
     *  This method is used in IMixtureStatModel::initializeStep
     *  in order to give a value to IStatModelBase::nbFreeParameter_.
     *  @return the number of free parameters
     **/
    virtual int computeNbFreeParameters() const = 0;

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
    inline virtual void setParameters() {}
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
};

/* Constructor.
 * @param nbCluster,nbSample number of clusters and samples
 **/
inline IMixtureStatModel::IMixtureStatModel( int nbSample, int nbCluster)
                                   : IStatModelBase(nbSample)
                                   , nbCluster_(nbCluster)
                                   , pk_(nbCluster), tik_(nbSample, nbCluster)
                                   , nk_(nbCluster), zi_(nbSample)
{}

/* copy constructor */
inline IMixtureStatModel::IMixtureStatModel( IMixtureStatModel const& model)
                                   : IStatModelBase(model)
                                   , nbCluster_(model.nbCluster_)
                                   , pk_(model.pk_), tik_(model.tik_)
                                   , nk_(model.nk_), zi_(model.zi_)
{}
/* destructor */
inline IMixtureStatModel::~IMixtureStatModel() {}

inline Real IMixtureStatModel::computeLikelihood(int i) const
{ return std::exp(computeLnLikelihood(i));}

/* @return the computed likelihood of the i-th sample.
 *  @param i index of the sample
 **/
inline Real IMixtureStatModel::computeLnLikelihood(int i) const
{
  // get maximal value
  CPointX lnComp(pk_.size());
  for (int k = pk_.begin(); k< pk_.end(); ++k)
  { lnComp[k] = std::log(pk_[k]) + lnComponentProbability(i, k);}
  // compute result
  Real lnCompMax = lnComp.maxElt();
  return std::log((lnComp-lnCompMax).exp().sum())+lnCompMax;
}

/* @return the computed log-likelihood. */
inline Real IMixtureStatModel::computeLnLikelihood() const
{
  Real res = 0.0;
  for (int i = tik().beginRows(); i< tik().endRows(); ++i)
  { res += computeLnLikelihood(i);}
  return res;
}

/* @return the computed ICL criteria. */
inline Real IMixtureStatModel::computeICL() const
{
  Real res = 0.0;
  for (int j = tik().beginCols(); j< tik().endCols(); ++j)
  { res += (tik_.col(j) * (tik_.col(j)+1e-15).log()).sum();}
  // compute result
  return (- 2. * lnLikelihood() + nbFreeParameter() * lnNbSample() - 2. * res);
}

/* @brief Initialize the model before its first use.
 *  This function can be overloaded in derived class for initialization of
 *  the specific model parameters. It should be called prior to any used of
 *  the class.
 *  @sa IMixture,MixtureBridge,MixtureLearner
 **/
inline void IMixtureStatModel::initializeStep()
{ setLnLikelihood(-Arithmetic<Real>::infinity());}

} // namespace STK

#endif /* STK_IMIXTURESTATMODEL_H */


