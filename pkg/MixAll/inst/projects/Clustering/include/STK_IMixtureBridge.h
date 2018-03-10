/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2016 Serge Iovleff

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  stkpp::Clustering
 * created on: 15 nov. 2013
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IMixtureBridge.h
 *  @brief In this file we define the interface for all the bridge classes
 *  between the mixtures and the composer.
 **/

#ifndef STK_IMIXTUREBRIDGE_H
#define STK_IMIXTUREBRIDGE_H

#include "STK_IMixture.h"

#include <DManager/include/STK_DataBridge.h>

namespace STK
{

/** @ingroup Clustering
 *  @brief Interface base class for the bridges of the STK++ mixture.
 *
 *  This interface produce a default implementation to all the virtual
 *  functions required by the STK::IMixture interface by delegating
 *  to the bridged mixture or parametersHandler the computations.
 *
 *  The virtual methods to implement in derived class are
 *  @code
 *    virtual Derived* create() const;
 *    virtual Derived* clone() const;
 *  @endcode
 *
 *  The members of this interface are the following
 *  @code
 *    Mixture mixture_;           // concrete class deriving from STK::IMixtureDensity
 *    ParamHandler paramHandler_; // specialization of the STK::ParametersHandler struct
 *    DataBridge<Data>* p_data_;
 *  @endcode
 *
 */
template<class Derived>
class IMixtureBridge: public IMixture, public IRecursiveTemplate<Derived>
{
  public:
    // iterator on the index pair of missing values
    typedef std::vector< std::pair<int,int> >::const_iterator ConstIterator;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Data Data;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Type Type;
    typedef typename hidden::MixtureBridgeTraits<Derived>::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits<Derived>::ParamHandler ParamHandler;
    enum
    {
      // class of mixture
      idMixtureClass_ = hidden::MixtureBridgeTraits<Derived>::idMixtureClass_
    };

  protected:
    /** default constructor.
     *  @param p_data pointer on the DataBridge wrapping the data set
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    IMixtureBridge( DataBridge<Data>* p_data, std::string const& idData, int nbCluster);
    /** copy constructor
     *  @param bridge the IMixtureBridge to copy
     **/
    IMixtureBridge( IMixtureBridge const& bridge);
    /** destructor */
    virtual ~IMixtureBridge() {}

  public:
    /** @return the mixture */
    inline Mixture const& mixture() const { return mixture_;}
    /** @return the structure managing the parameters. */
    inline ParamHandler const& paramHandler() const { return paramHandler_;}

    // getter and setter of the parameters
    /** set the parameters of the model */
    inline void setParameters( Parameters const& param) { mixture_.param_ = param;}
    /** get the parameters of the model */
    inline void getParameters( Parameters& param) const { param = mixture_.param_;}
    /** This function is used in order to set the current values of the
     *  parameters to the parameters struct of the mixture.
     *  @param param the array/expression with the parameters of the mixture to
     *  store in the ParamHandler.
     **/
    template<class Array>
    inline void setParameters(ExprBase<Array> const& param)
    { mixture_.param_.setParameters(param.asDerived());}
    /** This function is used in order to get the current values of the parameters
     *  in an array.
     *  @param param array that will store the parameters of the mixture.
     */
    template<class Array>
    inline void getParameters(Array& param) const { mixture_.getParameters(param);}

    // start default implementation of virtual method inherited from IMixture
    /** @brief Initialize the mixture model before its use by the composer.
     *  The parameters values are set to their default values if the mixture_
     *  is newly created. if IMixtureBridge::initializeStep is used during a
     *  cloning, mixture class have to take care of the existing values of the
     *  parameters.
     **/
    inline virtual void initializeStep()
    {
      if (!p_composer())
        STKRUNTIME_ERROR_NO_ARG(IMixtureBridge::initializeStep,composer is not set);
      if (!mixture_.initializeStep()) throw Clust::initializeStepFail_;
    }
     /** This function must be defined to return the component probability (PDF)
     *  for corresponding sample i and cluster k.
     *  @param i,k Sample and Cluster numbers
     *  @return the log-component probability
     **/
    inline virtual Real lnComponentProbability(int i, int k)
    { return mixture_.lnComponentProbability(i, k);}
    /** This function is equivalent to Mstep and must be defined to update
     *  parameters.
     **/
    inline virtual void paramUpdateStep()
    { if (!mixture_.run( p_tik(), p_tk())) throw Clust::mStepFail_;}
    /** @brief This function should be used in order to initialize randomly the
     *  parameters of the mixture.
     **/
    inline virtual void randomInit() { mixture_.randomInit( p_tik(), p_tk());};
    /** This function must return the number of free parameters.
     *  @return Number of free parameters
     **/
    inline virtual int nbFreeParameter() const { return mixture_.computeNbFreeParameters();}
    /** This function must return the number of variables.
     *  @return Number of variables
     **/
    inline virtual int nbVariable() const { return mixture_.nbVariable();}
    /** @brief This function should be used to store any intermediate results
     *  during various iterations after the burn-in period.
     *  @param iteration Provides the iteration number beginning after the burn-in
     *  period.
     **/
    inline virtual void storeIntermediateResults(int iteration)
    { paramHandler_.updateStatistics(mixture_.param_);}
    /**@brief This step can be used to signal to the mixtures that they must
     * release the stored results. This is usually called if the estimation
     * process failed.
     **/
    inline virtual void releaseIntermediateResults()
    { paramHandler_.releaseStatistics();}
    /** @brief set the parameters of the model.
     *  This step should be used to set the initialize the statistics to the
     *  current value of the parameters.
     **/
    inline virtual void setParametersStep()
    { paramHandler_.setStatistics(mixture_.param_);}
    /** @brief This step can be used by developer to finalize any thing. It will
     *  be called only once after we finish running the estimation algorithm.
     */
    inline virtual void finalizeStep() { mixture_.finalizeStep();}
    /** @brief This function should be used for Imputation of data.
     *  The default implementation (in the base class) is to do nothing.
     */
    inline virtual void imputationStep();
    /** @brief This function must be defined for simulation of all the latent
     * variables and/or missing data excluding class labels. The class labels
     * will be simulated by the framework itself because to do so we have to
     * take into account all the mixture laws.
     */
    inline virtual void samplingStep();
    /** This function can be used to write summary of parameters to the output stream.
     *  @param os Stream where you want to write the summary of parameters.
     */
    inline virtual void writeParameters(ostream& os) const
    { mixture_.writeParameters(p_tik(), os);}

  protected:
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    IMixtureBridge( Mixture const& mixture, std::string const& idData, int nbCluster);
    /** This function will be used for the imputation of the missing data
     *  at the initialization (initializationStep). By defalt missing values
     *  in a column will be replaced by the mean or map (maximum a
     *  posteriori) value of the column.
     **/
    void removeMissing();
    /** The Mixture to bridge with the composer */
    Mixture mixture_;
    /** The parameter handler associated to the mixture_ */
    ParamHandler paramHandler_;
    /** pointer on the data manager */
    DataBridge<Data>* p_data_;

};

/* copy constructor */
template< class Derived>
IMixtureBridge<Derived>::IMixtureBridge( IMixtureBridge const& bridge)
                                       : IMixture(bridge)
                                       , mixture_(bridge.mixture_)
                                       , paramHandler_(bridge.paramHandler_)
                                       , p_data_(bridge.p_data_)
{}
/* default constructor.
 *  @param p_data pointer on the DataBridge that will be used by the bridge.
 *  @param idData id name of the mixture model
 *  @param nbCluster number of cluster
 **/
template< class Derived>
IMixtureBridge<Derived>::IMixtureBridge( DataBridge<Data>* p_data, std::string const& idData, int nbCluster)
                                       : IMixture(idData)
                                       , mixture_(nbCluster)
                                       , paramHandler_(nbCluster)
                                       , p_data_(p_data)
{}
// implementation
template< class Derived>
void IMixtureBridge<Derived>::imputationStep()
{
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  { p_data_->elt(it->first, it->second) = mixture_.impute(it->first, it->second, p_tik()->row(it->first) );}
}
// implementation
template< class Derived>
void IMixtureBridge<Derived>::samplingStep()
{
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  { p_data_->elt(it->first, it->second) = mixture_.sample(it->first, it->second, p_tik()->row(it->first));}
}

// implementation
/* protected constructor to use in order to create a bridge.
 *  @param mixture the mixture to copy
 *  @param idData id name of the mixture
 *  @param nbCluster number of cluster
 **/
template< class Derived>
IMixtureBridge<Derived>::IMixtureBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                                       : IMixture(idData)
                                       , mixture_(mixture)
                                       , paramHandler_(nbCluster)
                                       , p_data_(0)
{}

template< class Derived>
void IMixtureBridge<Derived>::removeMissing()
{
  Type value = Type();
  int j, old_j = Arithmetic<int>::NA();
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  {
    j = it->second; // get column
    if (j != old_j)
    {
      old_j = j;
      value =  this->asDerived().safeValue(j);
    }
    p_data_->dataij()(it->first, j) = value;
  }
}


} // namespace STK

#endif /* IMIXTUREBRIDGE_H */
