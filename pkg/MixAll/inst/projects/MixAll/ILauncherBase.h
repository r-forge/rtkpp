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
 * created on: 10 May 2016
 * Author:   Iovleff, serge.iovleff@stkpp.org
 **/

/** @file ILauncherBase.h
 *  @brief In this file we define the interface base class for launchers.
 **/


#ifndef STK_ILAUNCHERBASE_H
#define STK_ILAUNCHERBASE_H

#include "RDataHandler.h"

namespace STK
{

/** The ILauncherBase allow to create the composer or learner for estimate, learn
 *  or predict a mixture model with less effort.
 **/
class ILauncherBase: public IRunnerBase
{
  public:
    /** constructor with a list of component.
     *  @param model a reference on the current model
     **/
    ILauncherBase( Rcpp::S4 model);
    /** destructor. */
    virtual ~ILauncherBase();
    /** @return the model */
    inline Rcpp::S4 const& s4_model() const { return s4_model_;}

  protected:
    /** set the parameters of a component given by its ID
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     **/
    template<class Derived>
    void setParametersToComponent( IMixtureStatModel* p_model
                                 , IMixtureManager<Derived> const& manager
                                 , String const& idData
                                 , Rcpp::S4 s4_component
                                 );
    /** set the parameters to a component given by its ID
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     **/
    void setParametersToComponent( IMixtureStatModel* p_model
                                 , KernelMixtureManager const& manager
                                 , String const& idData
                                 , Rcpp::S4 s4_component
                                 );
    /** set diagonal Gaussian parameters and data to S4 component
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setDiagGaussianParametersToComponent( IMixtureStatModel* p_model
                                             , DiagGaussianMixtureManager<RDataHandler> const& manager
                                             , String const& idData
                                             , Rcpp::S4 s4_component
                                             );
    /** set Poisson parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setPoissonParametersToComponent( IMixtureStatModel* p_model
                                        , PoissonMixtureManager<RDataHandler> const& manager
                                        , String const& idData
                                        , Rcpp::S4 s4_component
                                        );
    /** set gamma parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     *  */
    void setGammaParametersToComponent( IMixtureStatModel* p_model
                                      , GammaMixtureManager<RDataHandler> const& manager
                                      , String const& idData
                                      , Rcpp::S4 s4_component
                                      );
    /** set categorical parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setCategoricalParametersToComponent( IMixtureStatModel* p_model
                                            , CategoricalMixtureManager<RDataHandler> const& manager
                                            , String const& idData
                                            , Rcpp::S4 s4_component
                                            );
    /** set kernel parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setKernelParametersToComponent( IMixtureStatModel* p_model
                                       , KernelMixtureManager const& manager
                                       , String const& idData
                                       , Rcpp::S4 s4_component
                                       );

    /** get parameters of a component given by its ID
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     **/
    ArrayXX getParameters( Rcpp::S4 s4_component
                         , String const& idData
                         );
    /** get diagonal Gaussian parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getDiagGaussianParameters( Rcpp::S4 s4_component
                                     , String const& idData
                                     );
    /** get Poisson parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getPoissonParameters( Rcpp::S4 s4_component
                                , String const& idData
                                );
    /** get gamma parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     *  */
    ArrayXX getGammaParameters( Rcpp::S4 s4_component
                              , String const& idData
                              );
    /** get categorical parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getCategoricalParameters( Rcpp::S4 s4_component
                                    , String const& idData
                                    );
    /** get kernel parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getKernelParameters( Rcpp::S4 s4_component
                               , String const& idData
                               );

    /** model from the R side */
    Rcpp::S4 s4_model_;

};

template<class Derived>
void ILauncherBase::setParametersToComponent( IMixtureStatModel* p_model
                                            , IMixtureManager<Derived> const& manager
                                            , String const& idData
                                            , Rcpp::S4 s4_component)
{
  String rModelName = s4_component.slot("modelName");
  bool freeProp;
  switch (Clust::mixtureToMixtureClass(Clust::stringToMixture(rModelName, freeProp)))
  {
    case Clust::DiagGaussian_:
      setDiagGaussianParametersToComponent(p_model, manager.asDerived(), idData, s4_component);
      break;
    case Clust::Poisson_:
      setPoissonParametersToComponent(p_model, manager.asDerived(), idData, s4_component);
      break;
    case Clust::Gamma_:
      setGammaParametersToComponent(p_model, manager.asDerived(), idData, s4_component);
      break;
    case Clust::Categorical_:
      setCategoricalParametersToComponent(p_model, manager.asDerived(), idData, s4_component);
      break;
    case Clust::Kmm_:
      setKernelParametersToComponent(p_model, manager.asDerived(), idData, s4_component);
      break;
    case Clust::unknown_mixture_class_:
      break;
    default:
      break;
  }
}


} // namespace STK

#endif /* STK_ILAUNCHERBASE_H */
