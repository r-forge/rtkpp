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
#include "../inst/projects/MixAll/MixAll_Util.h"
#include "../inst/projects/MixAll/ClusterFacade.h"

using namespace Rcpp;

namespace STK
{
/* facade design pattern.
 * The KmmLauncher allow to create the strategy for estimating a mixture model
 * with less effort
 **/
KmmLauncher::KmmLauncher( Rcpp::S4 s4_model
                        , STK::KernelHandler const& handler
                        , Rcpp::IntegerVector const& nbCluster
                        , Rcpp::CharacterVector const& models)
                        : IRunnerBase()
                        , kernelHandler_(handler)
                        , kernelManager_(kernelHandler_)
                        , s4_model_(s4_model)
                        , s4_strategy_(s4_model_.slot("strategy"))
                        , v_models_(models)
                        , v_nbCluster_(nbCluster)
                        , criterion_(Rcpp::as<std::string>(s4_model_.slot("criterionName")))
                        , isMixedData_(false)
                        , p_composer_(0)
{}
///* facade design pattern.
// * The KmmLauncher allow to create the strategy for estimating a mixture model
// * with less effort
// **/
//KmmLauncher::KmmLauncher( SEXP model, SEXP nbCluster, SEXP critName )
//                        : IRunnerBase()
//                        , handlerDouble_()
//                        , handlerInt_()
//                        , s4_model_(model)
//                        , v_models_()
//                        , s4_strategy_(s4_model_.slot("strategy"))
//                        , handlerDouble_()
//                        , handlerInt_()
//                        , kernelManager_(kernelHandler_)
//{}
/* destructor. */
KmmLauncher::~KmmLauncher()
{}


/* create the mixtures in the given learner */
void KmmLauncher::createMixtures(IMixtureStatModel* p_model)
{
  p_model->createMixture(kernelManager_);
}

/* create the kernel mixtures in the given composer */
void KmmLauncher::updateMixtures(MixtureComposer* p_composer)
{
  typedef MixtureComposer::ConstMixtIterator ConstMixtIterator;
  Rcpp::S4 s4_component = s4_model_.slot("component");
  // loop over mixtures
  for (ConstMixtIterator it =  p_composer->v_mixtures().begin(); it != p_composer->v_mixtures().end(); it++)
  {
    std::string idData = (*it)->idData();
    std::string idModel;
    Clust::Mixture typeModel = Clust::stringToMixture(idModel);
    RVector<double> dim((SEXP)s4_component.slot("dim"));
    double kdim = (dim.size()>0) ? dim[0] : 10;
    kernelManager_.setDim(p_composer->getMixture(idData), kdim);
  }
}



/* get the kernel parameters */
void KmmLauncher::getParameters(IMixtureStatModel* p_model, std::string const& idData, Rcpp::S4 s4_component)
{
  // get parameters
  ArrayXX param;
  p_model->getParameters(kernelManager_,idData, param);
  // save results in s4_model
   s4_component.slot("sigma") = Rcpp::wrap(param.col(0));
  s4_component.slot("dim")   = Rcpp::wrap(param.col(1));
  // get data -- not necessary for kernels--
  //s4_component.slot("data") = kernelManager_.getData<double>(idData).matrix();
}

/* run the estimation */
bool KmmLauncher::run()
{

  // compute the best model
  Real initCriter = s4_model_.slot("criterion");
  Real criter = (isMixedData_) ? selectBestMixedModel()
                               : selectBestSingleModel();

  // get result common part of the estimated model
  s4_model_.slot("criterion")      = criter;
  s4_model_.slot("nbCluster")      = p_composer_->nbCluster();
  s4_model_.slot("lnLikelihood")   = p_composer_->lnLikelihood();
  s4_model_.slot("nbFreeParameter")= p_composer_->nbFreeParameter();
  s4_model_.slot("pk")             = Rcpp::wrap(p_composer_->pk());
  s4_model_.slot("tik")            = Rcpp::wrap(p_composer_->tik());
  s4_model_.slot("zi")             = Rcpp::wrap(p_composer_->zi());
  NumericVector fi = s4_model_.slot("lnFi");
  NumericVector zi = s4_model_.slot("zi");
  for (int i=0; i< fi.length(); ++i)
  {
    fi[i] = p_composer_->computeLnLikelihood(i);
    zi[i] += (1 - baseIdx);  // set base 1 for the class labels
  }
  if (criter == initCriter || !Arithmetic<Real>::isFinite(criter)) return false;
  return true;
}

/* get the parameters */
Real KmmLauncher::selectBestSingleModel()
{
  std::string idDataBestModel;
  // component
  Rcpp::S4 s4_component = s4_model_.slot("component");

  double bestCritValue = s4_model_.slot("criterion");
  int nbSample     = s4_model_.slot("nbSample");
  IMixtureComposer*  p_current =0;
  IMixtureCriterion* p_criterion =0;

  // start computation
  try
  {
    // create criterion
    p_criterion = createCriterion(criterion_);
    // start the estimation process, should end with the best model according to
    // the criterion
    ClusterFacade facade(p_current);
    facade.createFullStrategy(s4_strategy_);
    // loop over all the models
    for (int l=0; l <v_models_.size(); ++l)
    {
      std::string idData = "model" + typeToString<int>(l);
      std::string idModel(as<std::string>(v_models_[l]));
      bool freeProp;
      Clust::Mixture model = Clust::stringToMixture(idModel, freeProp);
      // for the current model,
      for (int k=0; k <v_nbCluster_.length(); ++k)
      {
        int K = v_nbCluster_[k];
        // create composer
        if (freeProp) { p_current = new MixtureComposer(nbSample, K);}
        else          { p_current = new MixtureComposerFixedProp(nbSample, K);}

        // create current mixture and register it
        std::string idData = "model" + typeToString<int>(l);
        createMixtures(static_cast<MixtureComposer*>(p_current));

        // run estimation and get results if possible
        if (!facade.run()) { msg_error_ += facade.error();}
        // compute criterion and update model if necessary
        p_criterion->setModel(p_current);
        p_criterion->run();
        if (bestCritValue > p_criterion->value())
        {
          if (p_composer_) { std::swap(p_current, p_composer_);}
          else             { p_composer_ = p_current; p_current = 0;}
          s4_component.slot("modelName") = idModel;
          idDataBestModel = idData;
          bestCritValue = p_criterion->value();
        }
        // release current composer
        if (p_current) { delete p_current; p_current = 0;}
      }
    }
    // release
    delete p_criterion;
    // get specific parameters
    getParameters(p_composer_, idDataBestModel, s4_component);
    return bestCritValue;
  }
  catch (Exception const& e)
  {
    if (p_current) delete p_current;
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
  Rcpp::List s4_list =s4_model_.slot("ldata");
  Real criter =s4_model_.slot("criterion");
  int nbSample =s4_model_.slot("nbSample");
  // main pointer
  IMixtureComposer* p_current =0;
  IMixtureCriterion* p_criterion =0;
  try
  {
    // create facade and strategy
    ClusterFacade facade(p_current);
    facade.createFullStrategy(s4_strategy_);
    bool sameProp = true;
    // loop over the list of component and fill handler_
    for (int l=0; l <s4_list.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = s4_list[l];
      std::string idData  = "model" + typeToString<int>(l);
      std::string idModel = s4_component.slot("modelName");
      // register
      bool freeMixture;
      Clust::Mixture model = Clust::stringToMixture(idModel, freeMixture);
      Clust::MixtureClass classModel = Clust::mixtureToMixtureClass(model);
      // if one of the model is free proportion, then we use free proportion
      sameProp &= (!freeMixture);
    }
    // create criterion
    p_criterion = createCriterion(criterion_);
    // loop over the models
    for (int k=0; k <v_nbCluster_.length(); ++k)
    {
      int K = v_nbCluster_[k];
      // create composer
      if (!sameProp) { p_current = new MixtureComposer(nbSample, K);}
      else           { p_current = new MixtureComposerFixedProp(nbSample, K);}
      // create all mixtures
      createMixtures(static_cast<MixtureComposer*>(p_current));
      // run estimation and get results if possible
      if (!facade.run()) { msg_error_ += facade.error();}
      // compute criterion and update model if necessary
      p_criterion->setModel(p_current);
      p_criterion->run();
      if (criter > p_criterion->value())
      {
        if (p_composer_) { std::swap(p_current, p_composer_);}
        else             { p_composer_ = p_current; p_current = 0;}
        criter = p_criterion->value();
      }
      // release current composer
      if (p_current) { delete p_current; p_current = 0;}
    }
    // release
    delete p_criterion;
    // get parameters
    for (int l=0; l <s4_list.size(); ++l)
    {
      // component
      Rcpp::S4 s4_component = s4_list[l];
      // id of the data set and of the model
      std::string idData  = "model" + typeToString<int>(l);
      getParameters(p_composer_, idData, s4_component);
    }
    //
    return criter;
  }
  catch (Exception const& e)
  {
    if (p_current) delete p_current;
    if (p_criterion) delete p_criterion;
    ::Rf_error(e.error().c_str()) ;
  }
  // failed
  return Arithmetic<Real>::max();
}




} // namespace STK

