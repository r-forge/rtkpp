/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2013  Serge Iovleff

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
 * Project:  stkpp::
 * created on: 3 sept. 2013
 * Author:   iovleff, serge.iovleff@stkpp.org
 **/

/** @file STK_MixtureStrategy.cpp
 *  @brief In this file we implement the strategies for estimating mixture model.
 **/

#include "STKernel/include/STK_Exceptions.h"
#include "../include/STK_MixtureStrategy.h"
#include "../include/STK_MixtureInit.h"
#include "../include/STK_MixtureAlgo.h"
#include "../include/STK_IMixtureComposer.h"

namespace STK
{

/* copy constructor
 *  @param strategy the strategy to copy
 **/
IMixtureStrategy::IMixtureStrategy( IMixtureStrategy const& strategy)
                                  : IRunnerBase(strategy), nbTry_(strategy.nbTry_)
                                  , p_model_(strategy.p_model_)
                                  , p_init_(strategy.p_init_->clone())
{}

/* store a model in p_model_ if it is better.
 * @param p_otherModel the model to store
 **/
void IMixtureStrategy::storeModel(IMixtureComposer*& p_otherModel)
{ if (p_model_->lnLikelihood()<p_otherModel->lnLikelihood())
  { std::swap(p_model_, p_otherModel);}
}

/* destructor */
IMixtureStrategy::~IMixtureStrategy() { if (p_init_) delete p_init_;}

/* destructor */
SimpleStrategyParam::~SimpleStrategyParam()
{ if (p_algo_) delete p_algo_;}

/* destructor */
XemStrategyParam::~XemStrategyParam()
{
  if (p_shortAlgo_) delete p_shortAlgo_;
  if (p_longAlgo_) delete p_longAlgo_;
}

/* destructor */
FullStrategyParam::~FullStrategyParam()
{
  if (p_shortAlgo_) delete p_shortAlgo_;
  if (p_longAlgo_) delete p_longAlgo_;
}

/* run the simple strategy */
bool SimpleStrategy::run()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("-----------------------------------------------\n");
  stk_cout << _T("Entering SimpleStrategy::run() with: ")
           << _T("nbTry_ = ") << this->nbTry_ << _T("\n");
#endif
  try
  {
    IMixtureComposer* p_currentModel = p_model_->create();
    // initialize and run algo. break if success
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
      // intialize current model
      p_init_->setModel(p_currentModel);
      if (p_init_->run())
      {
#ifdef STK_MIXTURE_VERBOSE
        stk_cout << _T("iTry =") << iTry << _T(". In SimpleStrategy::run(), Init step terminated. Current model:\n");
        p_currentModel->writeParameters(stk_cout);
        stk_cout << _T("\n\n");
#endif
        // run the estimation algorithm
        p_param_->p_algo_->setModel(p_currentModel);
        if (p_param_->p_algo_->run())
        {
#ifdef STK_MIXTURE_VERBOSE
            stk_cout << _T("iTry =") << iTry << ". In SimpleStrategy::run(), Long run terminated. current model\n";
            p_currentModel->writeParameters(stk_cout);
            stk_cout << _T("\n\n");
#endif
          // Check if we get a better result and break
          storeModel(p_currentModel); break;
        }  // algo step
      } // init step
    } // iTry
    delete p_currentModel;
  } catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "SimpleStrategy::run() terminated. current model\n";
  p_model_->writeParameters(stk_cout);
  stk_cout << "-----------------------------------------------\n";
#endif
  if (p_model_->lnLikelihood() == -STK::Arithmetic<Real>::max())
  {
    msg_error_ = STKERROR_NO_ARG(In SimpleStrategy::run,All trials failed);
    return false;
  }
  return true;
}

/* run the xem strategy */
bool XemStrategy::run()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering FullStrategy::run() with:\n")
           << _T("nbTry_ = ") << this->nbTry_ << _T("\n")
           << _T("nbShortRun_ = ") << p_param_->nbShortRun_ << _T("\n");
#endif
  // initialize bestModel and bestLikelihood
  try
  {
    // the current model is used in the short runs
    IMixtureComposer* p_currentModel     = p_model_->create();
    IMixtureComposer* p_currentBestModel = p_model_->create();
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("try number = ") << iTry << _T("\n");
#endif
      // find best of the shortModel and save it in p_currentBestModel
      for (int iShortRun = 0; iShortRun < p_param_->nbShortRun_; ++iShortRun)
      {
        // initialize current model
        p_init_->setModel(p_currentModel);
        if (p_init_->run())
        {
          // perform short run on the current model
          p_param_->p_shortAlgo_->setModel(p_currentModel);
          p_param_->p_shortAlgo_->run();
          // if we get a better result, swap it with currentBestModel
          if( p_currentBestModel->lnLikelihood()<p_currentModel->lnLikelihood())
          { std::swap(p_currentModel, p_currentBestModel);}
        } // initialization
      } // iShortRun
      // in case nbShortRun_==0
      // try to initialize bestCurrentModel, otherwise go to a next try
      if (p_param_->nbShortRun_ == 0)
      {
        // initialize current model
        p_init_->setModel(p_currentBestModel);
        if (!p_init_->run())
        { continue; }// model not initialized, we go to the next trial
      }
#ifdef STK_MIXTURE_VERBOSE
      stk_cout << _T("iTry =") << iTry
               << _T(". In FullStrategy::run(), short run terminated. best model:\n");
      p_currentBestModel->writeParameters(stk_cout);
      stk_cout << _T("\n\n");
#endif
      // start a long run with the better model
      p_param_->p_longAlgo_->setModel(p_currentBestModel);
      if (p_param_->p_longAlgo_->run()) { storeModel(p_currentBestModel); break;}
#ifdef STK_MIXTURE_VERBOSE
            stk_cout << "In FullStrategy::run(), Long run terminated. current model\n";
            p_model_->writeParameters(stk_cout);
            stk_cout << _T("\n\n");
#endif
    } // end iTry
    delete p_currentBestModel;
    delete p_currentModel;
  } catch (Exception const& e)
  {
    msg_error_ = e.error();
    return false;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << "FullStrategy::run() terminated. current model\n";
  p_model_->writeParameters(stk_cout);
  stk_cout << "-----------------------------------------------\n";
#endif
  if (p_model_->lnLikelihood() == -STK::Arithmetic<Real>::max())
  {
    msg_error_ = STKERROR_NO_ARG(In FullStrategy::run,All trials failed);
    return false;
  }
  return true;
}

/* run the xem strategy */
bool FullStrategy::run()
{
  IMixtureComposer* p_currentModel     = 0;
  IMixtureComposer* p_currentBestModel = 0;
  IMixtureComposer* p_bestShortModel   = 0;
  try
  {
    // the current best model store the best result throughout the strategy
    p_currentModel     = p_model_->create();
    p_currentBestModel = p_model_->create();
    p_bestShortModel   = p_model_->create();
    //
    // Main loop. If the Full strategy success in estimating a model, the
    // iterations are stopped and the best model find is stored in p_model_
    for (int iTry = 0; iTry < nbTry_; ++iTry)
    {
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("iTry = ") << iTry << _T("\n");
#endif
      // initialize the model with default values
      p_bestShortModel->initializeStep();
      // in case nbShortRun_==0: initialize directly p_bestShortModel
      if (p_param_->nbShortRun_ == 0)
      { initStep(p_currentModel, p_bestShortModel);}
      else
      {
        for (int iShort=0; iShort < p_param_->nbShortRun_; ++iShort)
        {
          // (re)-initialize current best model
          p_currentBestModel->initializeStep();
          // perform nbInitRun_ initialization step and get the best result
          // in p_currentBestModel
          initStep(p_currentModel, p_currentBestModel);
          // perform short run on the current best model
          p_param_->p_shortAlgo_->setModel(p_currentBestModel);
          if (p_param_->p_shortAlgo_->run())
          {
            // if we get a better result, store it in p_bestShortModel
            if( p_bestShortModel->lnLikelihood()<p_currentBestModel->lnLikelihood())
            { std::swap(p_bestShortModel, p_currentBestModel);}
          }
#ifdef STK_MIXTURE_VERBOSE
          stk_cout << _T("-------------------------------\n")
                   << _T("  iShort = ") << iShort << _T("\n")
                   << _T("p_bestShortModel->lnLikelihood() = ")
                   << p_bestShortModel->lnLikelihood() << _T("\n");
#endif
        } // ishort
      }
      // start a long run with p_bestShortModel. If success, save model
      // and exit the iTry loop
      p_param_->p_longAlgo_->setModel(p_bestShortModel);
      if (p_param_->p_longAlgo_->run())
      {
        // if we get a better result, store it in p_bestShortModel
        if( p_model_->lnLikelihood()<p_bestShortModel->lnLikelihood())
        { std::swap(p_model_, p_bestShortModel); break;}
      }
    } // end iTry
    // release memory
    delete p_currentBestModel; p_currentBestModel=0;
    delete p_bestShortModel; p_bestShortModel =0;
    delete p_currentModel; p_currentModel =0;
  }
  catch (Exception const& e)
  {
    if (p_currentBestModel) delete p_currentBestModel;
    if (p_bestShortModel)   delete p_bestShortModel;
    if (p_currentModel)     delete p_currentModel;
    msg_error_ = e.error();
    return false;
  }
  if (p_model_->lnLikelihood() == -STK::Arithmetic<Real>::max())
  {
    msg_error_ = STKERROR_NO_ARG(In FullStrategy::run,All trials failed);
    return false;
  }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("FullStrategy::run(), terminated without error.\n")
           << _T("Results:\n");
  p_model_->writeParameters(std::cout);
  stk_cout << _T("-------------------------------\n");
#endif
  return true;
}

/* Perform the Initialization step*/
void FullStrategy::initStep(IMixtureComposer*& p_currentModel, IMixtureComposer*& p_currentBestModel)
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("-------------------------------\n")
           << _T("Entering initStep\n");
#endif
  // set current model
  p_init_->setModel(p_currentModel);
  //
  // perform nbInitRun_ initialization (should be > 0)
  // and select the best model
  for (int iInitRun=0; iInitRun < p_param_->nbInitRun_; iInitRun++)
  {
    if (p_init_->run()) // Save only if we have a success
    {
      // if we get a better result, swap it with currentBestModel
      if( p_currentBestModel->lnLikelihood()<p_currentModel->lnLikelihood())
      { std::swap(p_currentModel, p_currentBestModel);}
    }
#ifdef STK_MIXTURE_VERBOSE
  stk_cout<< _T("     iInitRun = ") << iInitRun << _T("\n")
          << _T("p_currentBestModel->lnLikelihood() = ")
          << p_currentBestModel->lnLikelihood() << _T("\n");
#endif
  } // iInitRun
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Exiting initStep\n")
           << _T("-------------------------------\n");
#endif
}

} // namespace STK



