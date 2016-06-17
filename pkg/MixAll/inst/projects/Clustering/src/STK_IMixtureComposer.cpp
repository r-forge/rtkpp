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

/** @file STK_IMixtureComposer.h
 *  @brief In this file we implement the abstract base class for mixture models.
 **/

#include <cmath>

#ifdef STK_MIXTURE_DEBUG
#include <Arrays/include/STK_Display.h>
#endif

#include "../include/STK_IMixtureComposer.h"

#include <STatistiK/include/STK_Law_Categorical.h>
#include <STatistiK/include/STK_Stat_Functors.h>

namespace STK
{
IMixtureComposer::IMixtureComposer( int nbSample, int nbCluster)
                                  : IMixtureStatModel(nbSample, nbCluster)
                                  , state_(Clust::modelCreated_)
#ifndef _OPENMP
                                  , lnComp_(nbCluster)
#endif
{}

/* copy constructor */
IMixtureComposer::IMixtureComposer( IMixtureComposer const& model)
                                  : IMixtureStatModel(model)
                                  , state_(model.state_)
#ifndef _OPENMP
                                  , lnComp_(model.lnComp_)
#endif
{}
/* destructor */
IMixtureComposer::~IMixtureComposer() {}

/* @brief Initialize the model before at its first use.
 *  This function can be overloaded in derived class for initialization of
 *  the specific model parameters. It should be called prior to any used of
 *  the class.
 *  @sa IMixture,MixtureBridge,MixtureComposer
 **/
void IMixtureComposer::initializeStep()
{
  // (re)initialize the mixture parameters tik and pk. (virtual method)
  initializeMixtureParameters();
  // compute nk
  nk_  = Stat::sumByCol(tik_);
  // (re) initialize mixtures
  IMixtureStatModel::initializeStep();
  // compute zi (virtual method)
  mapStep();
  // compute proportions pk_ (virtual method)
  pStep();
  // (re)initialize the likelihood
  setLnLikelihood(computeLnLikelihood());
  // compute the number of free parameters
  setNbFreeParameter(computeNbFreeParameters());
  // compute number of free parameters and of variables
  setNbVariable(computeNbVariables());
  // update state
  setState(Clust::modelInitialized_);
}

/* initialize randomly the labels zi of the model */
void IMixtureComposer::randomClassInit()
{
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("Entering IMixtureComposer::randomClassInit()\n");
#endif
  // initialize mixture model if necessary
  if (state() < Clust::modelInitialized_) { initializeStep();}
  // generate random zi, compute tik and nk
  if (randomZi()<2) throw(Clust::randomClassInitFail_);
  // update parameters
  paramUpdateStep();
  // compute eStep()
  eStep();
  // model initialized
  setState(Clust::modelParamInitialized_);
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("IMixtureComposer::randomClassInit() done\n");
#endif
}

/* initialize randomly the posterior probabilities tik of the model */
void IMixtureComposer::randomFuzzyInit()
{
#ifdef STK_MIXTURE_VERBOSE
  stk_cout << _T("Entering IMixtureComposer::randomFuzzyInit(). state= ") << state() << _T("\n");
#endif
  // initialize mixture model if necessary
  if (state() < Clust::modelInitialized_) { initializeStep();}
  // create random tik and compute nk
  if (randomTik()<2) throw(Clust::randomFuzzyInitFail_);
  // compute zi
  mapStep();
  // update paramters values
  paramUpdateStep();
  // eStep
  eStep();
  // model intialized
  setState(Clust::modelParamInitialized_);
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("IMixtureComposer::randomFuzzyInit() done\n");
#endif
}

/* cStep */
int IMixtureComposer::cStep()
{
  for (int i=tik_.beginRows(); i < tik_.endRows(); i++)
  { cStep(i);}
  // count the number of individuals in each class
  nk_= Stat::sum(tik_);
  return nk_.minElt();
}

/* simulate zi  */
int IMixtureComposer::sStep()
{
  // simulate zi
  for (int i = zi_.begin(); i< zi_.end(); ++i)
  { sStep(i);}
  return cStep();
}

/* compute tik, default implementation. */
Real IMixtureComposer::eStep()
{
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("Entering IMixtureComposer::eStep()\n");
#endif
  Real sum = 0.; nk_ =0.;
  int i;
#ifdef _OPENMP
#pragma omp parallel for reduction (+:sum)
#endif
  for (i = tik_.beginRows(); i < tik_.endRows(); ++i) { sum += eStep(i);}
  // update ln-likelihood
  setLnLikelihood(sum);
  // compute proportions
  nk_ = Stat::sum(tik_);
#ifdef STK_MIXTURE_DEBUG
  stk_cout << _T("IMixtureComposer::eStep() done\n");
  stk_cout << _T("lnLikelihood =") << sum << _T("\n");
#endif
  return nk_.minElt();
}

/* Compute Zi using the Map estimate, default implementation. */
void IMixtureComposer::mapStep()
{
  for (int i = zi_.begin(); i< zi_.end(); ++i)
  { mapStep(i);  }
}

/* Simulate zi accordingly to tik.
 *  @param i index of the the individual
 **/
void IMixtureComposer::sStep(int i)
{ zi_.elt(i) = Law::Categorical::rand(tik_.row(i));}
/* Replace tik by zik
 *  @param i index of the the individual
 **/
void IMixtureComposer::cStep(int i)
{ tik_.row(i) = 0.; tik_.elt(i, zi_[i]) = 1.;}

/* compute tik, default implementation. */
Real IMixtureComposer::eStep(int i)
{
#ifdef _OPENMP
  CPointX lnComp_(tik_.cols());
#endif
  // get maximal element of ln(x_i,\theta_k) + ln(p_k)
  for (int k=tik_.beginCols(); k< tik_.endCols(); k++)
  { lnComp_[k] = std::log(pk_[k])+lnComponentProbability(i,k);}
  int kmax;
  Real max = lnComp_.maxElt(kmax);
  // set zi_
  zi_[i] = kmax;
  // return  max + sum_k p_k exp{lnCom_k - lnComp_kmax}
  Real sum =  (tik_.row(i) = (lnComp_ - max).exp()).sum();
  tik_.row(i) /= sum;
  return max + std::log( sum );
}

/* Compute Zi using the Map estimate, default implementation. */
void IMixtureComposer::mapStep(int i)
{
  int k;
  tik_.row(i).maxElt(k);
  zi_[i] = k;
}

/* estimate the proportions and the parameters of the components of the
 *  model given the current tik/zi mixture parameters values.
 **/
void IMixtureComposer::paramUpdateStep()
{ pStep();
  /* implement specific parameters estimation in concrete class. */
}

/* Compute prop using the ML estimate, default implementation. */
void IMixtureComposer::pStep()
{ pk_ = Stat::meanByCol(tik_);}

/* @brief Finalize the estimation of the model.
 *  The default behavior is compute current lnLikelihood.
 **/
void IMixtureComposer::finalizeStep()
{
  setLnLikelihood(computeLnLikelihood());
  setState(Clust::modelFinalized_);
}


// protected methods----------------------
/* Create the parameters of the  mixture model. */
void IMixtureComposer::initializeMixtureParameters()
{
  pk_  = 1./nbCluster_;
  tik_ = 1./nbCluster_;
}

/* generate random tik_ */
int IMixtureComposer::randomTik()
{
  nk_ = 0.;
  tik_.randUnif();
  for (int i = tik_.beginRows(); i < tik_.endRows(); ++i)
  {
    // create a reference on the i-th row
    CPointX tikRowi(tik_.row(i), true);
    tikRowi = tikRowi * pk_;
    tikRowi /= tikRowi.sum();
    nk_ += tikRowi;
  }
  return nk_.minElt();
}

/* generate random tik_ */
int IMixtureComposer::randomZi()
{
  Law::Categorical law(pk_);
  zi_.rand(law);
  return cStep();
}


} // namespace STK

