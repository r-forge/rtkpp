/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff, University Lille 1, Inria

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

    Contact : S..._DOT_I..._AT_stkpp.org (see copyright for ...)
*/

/*
 * Project:  MixAll
 * created on: 15 may 2016
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file KmmLauncher.cpp
 *  @brief In this file we implement the KmmLauncher which
 *  construct properly a mixture model.
 **/


#include "../inst/projects/MixAll/KmmLauncher.h"
#include "../inst/projects/MixAll/ClusterFacade.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The KmmLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
KmmLauncher::KmmLauncher( Rcpp::S4 s4_model
                        , Rcpp::IntegerVector const& nbCluster
                        , Rcpp::CharacterVector const& models

                        )
                        : ILauncherBase(s4_model)
                        , v_models_(models)
                        , v_nbCluster_(nbCluster)
                        , s4_strategy_(s4_model_.slot("strategy"))
                        , criterion_(Rcpp::as<std::string>(s4_model_.slot("criterionName")))
                        , isMixedData_(false)
                        , facade_()
{}
/* facade design pattern.
 * The KmmLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
KmmLauncher::KmmLauncher( Rcpp::S4 s4_model, Rcpp::IntegerVector const& nbCluster)
                        : ILauncherBase(s4_model)
                        , v_models_()
                        , v_nbCluster_(nbCluster)
                        , s4_strategy_(s4_model_.slot("strategy"))
                        , criterion_(Rcpp::as<std::string>(s4_model_.slot("criterionName")))
                        , isMixedData_(true)
                        , facade_()
{}
/* destructor. */
KmmLauncher::~KmmLauncher()
{}



/* run the estimation */
bool KmmLauncher::run()
{

  // compute the best model
  Real initCriter = s4_model_.slot("criterion");
  Real criter = (isMixedData_) ? selectBestMixedModel() : selectBestSingleModel();
  if (!facade_.p_composer_) return false;
  // get result common part of the estimated model
  s4_model_.slot("criterion")      = criter;
  s4_model_.slot("nbCluster")      = facade_.p_composer_->nbCluster();
  s4_model_.slot("lnLikelihood")   = facade_.p_composer_->lnLikelihood();
  s4_model_.slot("nbFreeParameter")= facade_.p_composer_->nbFreeParameter();
  s4_model_.slot("pk")             = Rcpp::wrap(facade_.p_composer_->pk());
  s4_model_.slot("tik")            = Rcpp::wrap(facade_.p_composer_->tik());
  s4_model_.slot("zi")             = Rcpp::wrap(facade_.p_composer_->zi());
  NumericVector fi = s4_model_.slot("lnFi");
  IntegerVector zi = s4_model_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = facade_.p_composer_->computeLnLikelihood(i);
    zi[i] += (1 - baseIdx);  // set base 1 for the class labels
  }
  if (criter == initCriter || !Arithmetic<Real>::isFinite(criter)) return false;
  return true;
}


/* Create kernels and add them to KernelHandler from facade
 **/
Kernel::IKernel* KmmLauncher::createKernel(Rcpp::S4 s4_component, std::string const& idData, KernelMixtureManager& manager)
{
  // get kernel parameters from component
  std::string         kernelName        = Rcpp::as<std::string>(s4_component.slot("kernelName"));
  Rcpp::NumericVector kernelParameters  = s4_component.slot("kernelParameters");
  bool                kernelComputation = Rcpp::as<bool>(s4_component.slot("kernelComputation"));

  Kernel::IKernel* p_kernel = 0;
  // build kernel
  switch (Kernel::stringToKernelType(kernelName))
  {
    case Kernel::exponential_:
    {
      Real param1 = 1;
      if (kernelParameters.length()>0) { param1 = kernelParameters[0];}
      NumericMatrix r_data = s4_component.slot("data");
      DataBridge< RMatrix<double> >* p_bridge = new DataBridge< RMatrix<double> >( idData, RMatrix<double>(r_data));
      manager.registerDataBridge(p_bridge);
      p_kernel = new Kernel::Exponential< RMatrix<double> >(p_bridge->dataij(), param1);
      break;
    }
    case Kernel::gaussian_:
    {
      Real param1 = 1;
      if (kernelParameters.length()>0) { param1 = kernelParameters[0];}
      NumericMatrix r_data = s4_component.slot("data");
      DataBridge< RMatrix<double> >* p_bridge = new DataBridge< RMatrix<double> >( idData, RMatrix<double>(r_data));
      manager.registerDataBridge(p_bridge);
      p_kernel = new Kernel::Gaussian< RMatrix<double> >(p_bridge->dataij(), param1);
      break;
    }
    case Kernel::laplace_:
    {
      Real param1 = 1;
      if (kernelParameters.length()>0) { param1 = kernelParameters[0];}
      NumericMatrix r_data = s4_component.slot("data");
      DataBridge< RMatrix<double> >* p_bridge = new DataBridge< RMatrix<double> >( idData, RMatrix<double>(r_data));
      manager.registerDataBridge(p_bridge);
      p_kernel = new Kernel::Laplace< RMatrix<double> >(p_bridge->dataij(), param1);
      break;
    }
    case Kernel::linear_:
    {
      NumericMatrix r_data = s4_component.slot("data");
      DataBridge< RMatrix<double> >* p_bridge = new DataBridge< RMatrix<double> >( idData, RMatrix<double>(r_data));
      manager.registerDataBridge(p_bridge);
      p_kernel = new Kernel::Linear< RMatrix<double> >(p_bridge->dataij());
      break;
    }
    case Kernel::polynomial_:
    {
      Real param1 = 2, param2 =0;
      if (kernelParameters.length()>0) { param1 = kernelParameters[0];}
      if (kernelParameters.length()>2) { param2 = kernelParameters[1];}
      NumericMatrix r_data = s4_component.slot("data");
      DataBridge< RMatrix<double> >* p_bridge = new DataBridge< RMatrix<double> >( idData, RMatrix<double>(r_data));
      manager.registerDataBridge(p_bridge);
      p_kernel = new Kernel::Polynomial< RMatrix<double> >(p_bridge->dataij(), param1, param2);
      break;
    }
    case Kernel::rationalQuadratic_:
    {
      Real param1 = 1;
      if (kernelParameters.length()>0) { param1 = kernelParameters[0];}
      NumericMatrix r_data = s4_component.slot("data");
      DataBridge< RMatrix<double> >* p_bridge = new DataBridge< RMatrix<double> >( idData, RMatrix<double>(r_data));
      manager.registerDataBridge(p_bridge);
      p_kernel = new Kernel::RationalQuadratic< RMatrix<double> >(p_bridge->dataij(), param1);
      break;
    }
    case Kernel::hamming_:
    {
      Real param1 = 2, param2 =0;
      if (kernelParameters.length()>0) { param1 = kernelParameters[0];}
      if (kernelParameters.length()>2) { param2 = kernelParameters[1];}
      IntegerMatrix r_data = s4_component.slot("data");
      DataBridge< RMatrix<int> >* p_bridge = new DataBridge< RMatrix<int> >( idData, RMatrix<int>(r_data));
      manager.registerDataBridge(p_bridge);
      p_kernel = new Kernel::Hamming< RMatrix<int> >(p_bridge->dataij(), param1);
      break;
    }
    default:
      return 0;
      break;
  }
  // compute Gram matrix
  if (kernelComputation)
  {
    p_kernel->run();
    s4_component.slot("gram") = STK::wrap(p_kernel->gram());
  }
  return p_kernel;
}

/* get the parameters */
Real KmmLauncher::selectBestSingleModel()
{
  // component
  Rcpp::S4 s4_component = s4_model_.slot("component");
  double critBestModel  = s4_model_.slot("criterion");
  int nbSample          = s4_model_.slot("nbSample");
  Rcpp::NumericVector r_dim = s4_component.slot("dim");
  double dim = r_dim[0];

  //
  std::string idBestModel;
  bool freeProp;
  IMixtureComposer*  p_bestModel =0;
  IMixtureCriterion* p_criterion =0;

  // start computation
  try
  {
    // create criterion
    p_criterion = Clust::createCriterion(criterion_);
    if (!p_criterion)
    {
      msg_error_ = STKERROR_1ARG(KmmLauncher::run,criterion_,Erro in criterion creation);
      return Arithmetic<Real>::max();
    }

    // start the estimation process, should end with the best model according to the criterion
    ClusterFacade facade(facade_.p_composer_);
    facade.createFullStrategy(s4_strategy_);

    // loop over all the models
    for (int l=0; l <v_models_.length(); ++l)
    {
      std::string idData  = "model" + typeToString<int>(l);
      std::string idModel = as<std::string>(v_models_[l]);

      // create kernel and register it
      Kernel::IKernel* p_kernel = createKernel(s4_component, idData, facade_.kmmManager_);
      if (!facade_.kerHandler_.addKernel(p_kernel, idData, idModel))
      { msg_error_ = STKERROR_1ARG(KmmLauncher::run,idData,Error in kernel creation);
        return Arithmetic<Real>::max();
      }

      // get freeProp
      Clust::stringToMixture(idModel, freeProp);

      // loop over all number of cluster
      for (int k=0; k <v_nbCluster_.length(); ++k)
      {
        int nbCluster = v_nbCluster_[k];
        // create composer
        facade_.p_composer_ = (freeProp) ? new MixtureComposer(nbSample, nbCluster)
                                         : new MixtureComposerFixedProp(nbSample, nbCluster);
        IMixture* p_mixture = facade_.kmmManager_.createMixture(idData, nbCluster);
        facade_.kmmManager_.setDim(p_mixture, dim);
        facade_.p_composer_->registerMixture(p_mixture);

        // run estimation and get results if possible
        if (!facade.run())
        { msg_error_ += facade.error();}

        // compute criterion and update model if necessary
        p_criterion->setModel(facade_.p_composer_);
        p_criterion->run();
        if (critBestModel > p_criterion->value())
        {
          idBestModel   = idData;
          critBestModel = p_criterion->value();
          std::swap(p_bestModel, facade_.p_composer_);
          s4_component.slot("modelName") = idModel;
        }
        // release current composer as best model is pointed by p_bestModel
        if (facade_.p_composer_)
        { delete facade_.p_composer_; facade_.p_composer_ = 0;}
      }
    }
    // save best model
    facade_.p_composer_ = p_bestModel;
    if (facade_.p_composer_) // in case every thing fail
    {
      setKernelParametersToComponent( facade_.p_composer_, facade_.kmmManager(), idBestModel
                                    , s4_component);
    }
    if (p_criterion) { delete p_criterion;}
    return critBestModel;
  }
  catch (Exception const& e)
  {
    if (p_bestModel) delete p_bestModel;
    if (p_criterion) delete p_criterion;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}

/* select best mixed data model */
Real KmmLauncher::selectBestMixedModel()
{
  // list of the component
  Rcpp::List list_component = s4_model_.slot("lcomponent");
  double critBestModel      = s4_model_.slot("criterion");
  int nbSample              = s4_model_.slot("nbSample");
  // auxiliaries pointers
  IMixtureComposer*  p_bestModel =0;
  IMixtureCriterion* p_criterion =0;
  try
  {
    // create criterion
    p_criterion = Clust::createCriterion(criterion_);
    if (!p_criterion)
    {
      msg_error_ = STKERROR_1ARG(KmmLauncher::run,criterion_,Error in criterion creation);
      return Arithmetic<Real>::max();
    }

    // start the estimation process, should end with the best model according to the criterion
    ClusterFacade facade(facade_.p_composer_);
    facade.createFullStrategy(s4_strategy_);

    // loop over the list of component and fill handler_
    bool sameProp = true;
    for (int l=0; l <list_component.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = list_component[l];
      std::string idData  = "model" + typeToString<int>(l);
      std::string idModel = s4_component.slot("modelName");

      // create kernel and register it
      Kernel::IKernel* p_kernel = createKernel(s4_component, idData, facade_.kmmManager_);
      if (!facade_.kerHandler_.addKernel(p_kernel, idData, idModel))
      { msg_error_ = STKERROR_1ARG(KmmLauncher::run,idData,Error in kernel creation);
        return Arithmetic<Real>::max();
      }

      // compute sameProp. If one of the model is free proportion, then we use free proportion
      bool freeMixture;
      Clust::Mixture model = Clust::stringToMixture(idModel, freeMixture);
      sameProp &= (!freeMixture);
    }

    // loop over the number of clusters
    for (int k=0; k <v_nbCluster_.length(); ++k)
    {
      int nbCluster = v_nbCluster_[k];
      // create composer
      if (!sameProp) { facade_.p_composer_ = new MixtureComposer(nbSample, nbCluster);}
      else           { facade_.p_composer_ = new MixtureComposerFixedProp(nbSample, nbCluster);}

      // create all mixtures
      for (int l=0; l <list_component.length(); ++l)
      {
        // get component and dimension to specific model
        Rcpp::S4 s4_component = list_component[l];
        Rcpp::NumericVector r_dim = s4_component.slot("dim");
        double dim = r_dim[0];
        std::string idData  = "model" + typeToString<int>(l);
        IMixture* p_mixture = facade_.kmmManager_.createMixture(idData, nbCluster);
        facade_.kmmManager_.setDim(p_mixture, dim);
        facade_.p_composer_->registerMixture(p_mixture);
      }

      // run estimation and get results if possible
      if (!facade.run()) { msg_error_ += facade.error();}
      // compute criterion and update model if necessary
      p_criterion->setModel(facade_.p_composer_);
      p_criterion->run();
      if (critBestModel > p_criterion->value())
      {
        critBestModel = p_criterion->value();
        std::swap(p_bestModel, facade_.p_composer_);
      }
      // release current composer as best model is pointed by p_bestModel
      if (facade_.p_composer_)
      { delete facade_.p_composer_; facade_.p_composer_ = 0;}
    }
    //
    if (p_criterion) { delete p_criterion;}
    //save best model
    if (p_bestModel) // in case everything failed
    {
      facade_.p_composer_ = p_bestModel;
      // get parameters
      for (int l=0; l <list_component.length(); ++l)
      {
        // component
        Rcpp::S4 s4_component = list_component[l];
        // id of the data set and of the model
        std::string idData  = "model" + typeToString<int>(l);
        // get specific parameters
        setKernelParametersToComponent( facade_.p_composer_
                                      , facade_.kmmManager_
                                      , idData
                                      , s4_component);
      }
    }
    return critBestModel;
  }
  catch (Exception const& e)
  {
    if (p_bestModel) delete p_bestModel;
    if (p_criterion) delete p_criterion;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}




} // namespace STK

