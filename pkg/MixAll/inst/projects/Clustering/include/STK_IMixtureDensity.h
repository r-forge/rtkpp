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
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_IMixtureDensity.h
 *  @brief In this file we define the main interface class for mixture models.
 **/


#ifndef STK_IMIXTUREDENSITY_H
#define STK_IMIXTUREDENSITY_H

#include "STK_Clust_Util.h"

#include <StatModels/include/STK_Model_Util.h>
#include <STatistiK/include/STK_Law_Categorical.h>

#ifdef STK_MIXTURE_DEBUG
#include "Arrays/include/STK_Display.h"
#endif

namespace STK
{
namespace hidden
{
/** Main class for the mixtures traits policy.
 *  The traits struct MixtureTraits must be specialized for any
 *  Mixture deriving from the Interface IMixtureDensity.
 **/
template <class Mixture> struct MixtureTraits;

} // namespace hidden

/** Parameters class. All statistical models has an unique Id defined
 *  in STK::Clust::Mixture enumeration.
 **/
template <int Id> struct ModelParameters;

/**@ingroup Clustering
 *  @brief Base class for all Mixture densities.
 *
 *  Let @e X be an arbitrary measurable space and let
 *  \f$ \mathbf{x} =\{{\mathbf x}_1,...,{\mathbf x}_n\}\f$
 *  be @e n independent vectors in @e X. If each \f${\mathbf x}_i\f$
 *  arises from a probability distribution with density
 *  \f[
 *    f({\mathbf x}_i|\theta) = \sum_{k=1}^K p_k h({\mathbf x}_{i}| \lambda_{k},\alpha)
 *  \f]
 *  where the \f$p_k\f$'s are the mixing proportions, \f$h(\cdot| \lambda_{k},\alpha)\f$
 *  denotes a distribution parameterized by \f$\lambda_k\f$ and \f$\alpha\f$, it
 *  is said that we observe a mixture model.
 *
 *  This interface class is the base class for all distributions part of a mixture
 *  model.
 *
 *  @sa IMixture, IMixtureDensity, MixtureComposer
 *
 * At this level the Parameters struct is created. The call to @e setData
 * trigger the call to @e initializeModel which itself trigger a call to
 * @e initializeModelImpl eventually implemented in the derived class.
 *
 * The pseudo virtual methods (needed by IMixtureBridge) to implement in derived
 * classes are
 * @code
 *   template<class Weights>
 *   Type impute(int i, int j, Weights const& pk) const;
 *   Type rand(int i, int j, int k) const;
 *   Real lnComponentProbability(int i, int k) const;
 *   void randomInit( CArrayXX const*  p_tik, CPointX const* p_nk) ;
 *   bool run( CArrayXX const*  p_tik, CPointX const* p_nk) ;
 *   int computeNbFreeParameters() const;
 * @endcode
 *
 * The pseudo virtual methods to implement if needed in derived classes are
 * @code
 *   // default implementation (do nothing) provided to all these methods
 *   void initializeModelImpl();
 *   bool initializeStepImpl(); // return true by default
 *   void finalizeStepImpl();
 * @endcode
 *
 * @sa IMixtureDensityBase, IRecursiveTemplate
 **/
template<class Derived>
class IMixtureDensity: public IRecursiveTemplate<Derived>
{
  public:
    typedef typename hidden::MixtureTraits<Derived>::Array Array;
    typedef typename hidden::MixtureTraits<Derived>::Parameters Parameters;
    typedef typename Array::Type Type;

  protected:
    /** Default constructor.
     *  @param nbCluster the number of cluster
     **/
    IMixtureDensity( int nbCluster);
    /** copy constructor.
     *  - The Parameter class is copied using the copy constructor.
     *  - The pointer on the data set is copied as-is. Check if you should not
     *  change it on the copied object.
     *  @param model the model to copy
     **/
    IMixtureDensity( IMixtureDensity const& model);

  public:
    /** destructor */
    inline ~IMixtureDensity() {}
    /** create pattern.  */
    inline IMixtureDensity* create() const { return new Derived(this->nbCluster());}

    /** @return the number of cluster */
    inline int nbCluster() const { return nbCluster_;}
    /** @return the total available observations */
    inline int nbSample() const { return nbSample_;}
    /** @return the Log of the total available observations */
    inline Real lnNbSample() const
    { return (nbSample_ <= 0) ? -Arithmetic<Real>::infinity() : std::log((Real)nbSample_);}
    /** @return the total available variables */
    inline int nbVariable() const { return nbVariable_;}
    /** @return a pointer on the current data set */
    inline Array const* p_data() const { return p_dataij_;}

    /** @brief Set the data set.
     *  Setting a (new) data set will trigger the initializeModel() method.
     *  @param data the data set to set
     **/
    void setData(Array const& data);
    /** @brief This function will be called at the beginning of the estimation
     *  process once the model is created and data is set.
     *  @note a stk++ mixture create and initialize all the containers when the data
     *  is set. Thus the default behavior is @c return true.
     */
    inline bool initializeStep() { return this->asDerived().initializeStepImpl();}
    /** set the parameters obtained with the intermediate results and release
     *  the intermediate results. */
    inline void setParametersStep()
    {
      param_.setParametersStep();
      this->asDerived().setParametersImpl();
    }
    /** @brief This function will be called once the model is estimated.
     *  perform specific model finalization stuff */
    inline void finalizeStep() { this->asDerived().finalizeStepImpl();}

    // default implementation of the pseudo-virtual methods
    /** default implementation of initializeModelImpl (do nothing) */
    inline void initializeModelImpl() {}
    /** default implementation of initializeStepImpl (return true) */
    inline bool initializeStepImpl() { return true;}
    /** default implementation of finalizeStepImpl (do nothing) */
    inline void finalizeStepImpl() {}

    /** @return a simulated value for the jth variable of the ith sample
     *  @param i,j indexes of the data to impute
     *  @param tk the probabilities of each class for the ith individual
     **/
    template<class Weights>
    inline Type sample(int i, int j, Weights const& tk) const
    { return this->asDerived().rand(i, j, Law::Categorical::rand(tk));}

  protected:
    /** @brief Initialize the model before its first use.
     * This function is triggered when data set is set.
     * In this interface, the @c initializeModel() method
     *  - set the number of samples and variables of the mixture model
     *  - call the derived class implemented method
     * @code
     *   initializeModelImpl()
     * @endcode
     * for initialization of the specific model parameters if needed.
     **/
    void initializeModel();

    /** Set the number of sample of the model
     *  @param nbSample number of sample of the model
     * */
    inline void setNbSample( int nbSample) { nbSample_ = nbSample;}
    /** Set the number of variables of the model
     *  @param nbVariable number of variables of the model
     * */
    inline void setNbVariable( int nbVariable) { nbVariable_ = nbVariable;}

    /** parameters of the derived mixture model */
    Parameters param_;

  private:
    /** number of cluster. */
    int nbCluster_;
    /** total available samples */
    int nbSample_;
    /** total available variables */
    int nbVariable_;
    /** pointer on the data set */
    Array const* p_dataij_;
};


/* Default constructor.
 *  @param nbCluster the number of cluster
 **/
template<class Derived>
IMixtureDensity<Derived>::IMixtureDensity( int nbCluster)
                                                : param_(nbCluster)
                                                , nbCluster_(nbCluster)
                                                , nbSample_(0)
                                                , nbVariable_(0)
                                                , p_dataij_(0)
{}
/* copy constructor.
 *  - The Parameter class is copied using the copy constructor.
 *  - The pointer on the data set is copied as-is. Check if you should not
 *  change it on the copied object.
 *  @param model the model to copy
 **/
template<class Derived>
IMixtureDensity<Derived>::IMixtureDensity( IMixtureDensity const& model)
             : param_(model.param_)
             , nbCluster_(model.nbCluster_)
             , nbSample_(model.nbSample_)
             , nbVariable_(model.nbVariable_)
             , p_dataij_(model.p_dataij_)
{}

/* @brief Set the data set.
 *  Setting a (new) data set will trigger the initializeModel() method.
 *  @param data the data set to set
 **/
template<class Derived>
void IMixtureDensity<Derived>::setData(Array const& data)
{
  p_dataij_ = &data;
  initializeModel();
}

/* @brief Initialize the model before its first use.
 * This function is triggered when data set is set.
 * In this interface, the @c initializeModel() method
 *  - set the number of samples and variables of the mixture model
 *  - call the derived class implemented method
 * @code
 *   initializeModelImpl()
 * @endcode
 * for initialization of the specific model parameters if needed.
 **/
template<class Derived>
void IMixtureDensity<Derived>::initializeModel()
{
  // set dimensions
  this->setNbSample(p_dataij_->sizeRows());
  this->setNbVariable(p_dataij_->sizeCols());
  // call specific model initialization stuff
  this->asDerived().initializeModelImpl();
}

} // namespace STK

#endif /* STK_IMIXTUREDENSITY_H */
