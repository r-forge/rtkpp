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
    typedef IMixtureManager<RDataHandler> IManager;
    /** constructor with a list of component.
     *  @param model a reference on the current model
     **/
    ILauncherBase( Rcpp::S4 model);
    /** destructor. */
    virtual ~ILauncherBase();
    /** @return the model */
    inline Rcpp::S4 const& s4_model() const { return s4_model_;}

  protected:
    /** get the parameters of a component given by its ID
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     **/
    void setParametersToComponent( IMixtureStatModel* p_model
                                 , IManager const& manager
                                 , std::string const& idData
                                 , Rcpp::S4 s4_component
                                 );
    /** get the parameters of a component given by its ID
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     **/
    void setParametersToComponent( IMixtureStatModel* p_model
                                 , KernelMixtureManager const& manager
                                 , std::string const& idData
                                 , Rcpp::S4 s4_component
                                 );
    /** set diagonal Gaussian parameters and data to S4 component
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setDiagGaussianParametersToComponent( IMixtureStatModel* p_model
                                             , IManager const& manager
                                             , std::string const& idData
                                             , Rcpp::S4 s4_component
                                             );
    /** get Poisson parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setPoissonParametersToComponent( IMixtureStatModel* p_model
                                        , IManager const& manager
                                        , std::string const& idData
                                        , Rcpp::S4 s4_component
                                        );
    /** get gamma parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     *  */
    void setGammaParametersToComponent( IMixtureStatModel* p_model
                                      , IManager const& manager
                                      , std::string const& idData
                                      , Rcpp::S4 s4_component
                                      );
    /** get categorical parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setCategoricalParametersToComponent( IMixtureStatModel* p_model
                                            , IManager const& manager
                                            , std::string const& idData
                                            , Rcpp::S4 s4_component
                                            );
    /** get kernel parameters
     *  @param p_model model with the parameters to get
     *  @param manager mixture manager
     *  @param idData ID of the data
     *  @param s4_component component storing the parameters
     * */
    void setKernelParametersToComponent( IMixtureStatModel* p_model
                                       , KernelMixtureManager const& manager
                                       , std::string const& idData
                                       , Rcpp::S4 s4_component
                                       );

    /** get parameters of a component given by its ID
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     **/
    ArrayXX getParameters( Rcpp::S4 s4_component
                         , std::string const& idData
                         );
    /** get diagonal Gaussian parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getDiagGaussianParameters( Rcpp::S4 s4_component
                                     , std::string const& idData
                                     );
    /** get Poisson parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getPoissonParameters( Rcpp::S4 s4_component
                                , std::string const& idData
                                );
    /** get gamma parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     *  */
    ArrayXX getGammaParameters( Rcpp::S4 s4_component
                              , std::string const& idData
                              );
    /** get categorical parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getCategoricalParameters( Rcpp::S4 s4_component
                                    , std::string const& idData
                                    );
    /** get kernel parameters
     *  @param s4_component the component with the parameters to get
     *  @param idData ID of the data
     * */
    ArrayXX getKernelParameters( Rcpp::S4 s4_component
                               , std::string const& idData
                               );

    /** model from the R side */
    Rcpp::S4 s4_model_;

};

} // namespace STK

#endif /* STK_ILAUNCHERBASE_H */
