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
 * Originally created by Parmeet Bhatia <b..._DOT_p..._AT_gmail_Dot_com>
 **/

/** @file STK_PoissonBridge.h
 *  @brief In this file we define the bridge class between the Poisson mixture
 *  models and the composer.
 **/

#ifndef STK_POISSONBRIDGE_H
#define STK_POISSONBRIDGE_H

#include "STK_MixturePoisson_ljk.h"
#include "STK_MixturePoisson_ljlk.h"
#include "STK_MixturePoisson_lk.h"

#include "../STK_IMixtureBridge.h"

#include <DManager/include/STK_DataBridge.h>
#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{

// forward declaration
template<int Id, class Data> class PoissonBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixturePoisson_ljk model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge<Clust::Poisson_ljk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef MixturePoisson_ljk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Poisson_ljk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Poisson_ljk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};

/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixturePoisson_lk model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge<Clust::Poisson_lk_, Data_> >
{
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef MixturePoisson_lk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Poisson_lk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Poisson_lk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixturePoisson_ljlk model
 **/
template<class Data_>
struct MixtureBridgeTraits< PoissonBridge< Clust::Poisson_ljlk_, Data_> >
{
  typedef Data_ Data;
  /** Type of the mixture model */
  typedef MixturePoisson_ljlk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Poisson_ljlk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Poisson_ljlk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Poisson_
  };
};

} // namespace hidden

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixturePoisson_ljlk model
 **/
template <>
struct ParametersHandler<Clust::Poisson_ljlk_>
{
    typedef ModelParameters<Clust::Poisson_ljlk_> Parameters;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<Real, Real> > stat_lambdak_;
    /** Array of the lambdaj_ statistics */
    Stat::Online<CVectorX, Real>  stat_lambdaj_;
    /** default constructor. All lambdas are initialized to 1. */
    inline ParametersHandler( int nbCluster)
                            : stat_lambdak_(nbCluster), stat_lambdaj_()
    {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_lambdak_(param.stat_lambdak_)
                            , stat_lambdaj_(param.stat_lambdaj_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_lambdaj_ = other.stat_lambdaj_;
      stat_lambdak_ = other.stat_lambdak_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
      { stat_lambdak_[k].update(param.lambdak_[k]);}
      stat_lambdaj_.update(param.lambdaj_);
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
      {
        param.lambdak_[k] = stat_lambdak_[k].mean();
        stat_lambdak_[k].release();
      }
      param.lambdaj_ = stat_lambdaj_.mean();
      stat_lambdaj_.release();

    }
    /** Get the parameters of the mixture model.
     *  It is assumed that array params stores for each class the lambdaj and
     *  lambdak parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.lambdak_.begin(); k<param.lambdak_.end(); ++k)
      {
        for(int j=param.lambdaj_.begin(); j<param.lambdaj_.end(); ++j)
        { params(k,j) = param.lambdaj_[j] * param.lambdak_[k];}
      }
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      CVectorX lambdak = Stat::meanByRow(params.asDerived());
      CPointX lambdaj  = Stat::meanByCol(params.asDerived());
      Real cte = std::sqrt((params/(lambdak * lambdaj)).mean());
      param.lambdak_ = cte * lambdak;
      param.lambdaj_ = cte * lambdaj;
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
      { stat_lambdak_[k].release();}
      stat_lambdaj_.release();
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      stat_lambdaj_.resize(range);
      for(int k=stat_lambdak_.begin(); k<stat_lambdak_.end(); ++k)
      { stat_lambdak_[k].release();}
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixturePoisson_ljk model
 **/
template <>
struct ParametersHandler<Clust::Poisson_ljk_>
{
    typedef ModelParameters<Clust::Poisson_ljk_> Parameters;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<CPointX, Real> > stat_lambda_;
    /** default constructor. All lambdas are initialized to 1. */
    inline ParametersHandler( int nbCluster)
                            : stat_lambda_(nbCluster)
    {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_lambda_(param.stat_lambda_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_lambda_ = other.stat_lambda_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      { stat_lambda_[k].update(param.lambda_[k]);}
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      {
        param.lambda_[k] = stat_lambda_[k].mean();
        stat_lambda_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.lambda_.begin(); k<param.lambda_.end(); ++k)
      { params.row(k) = param.lambda_;}
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      for(int k=param.lambda_.begin(); k<param.lambda_.end(); ++k)
      { param.lambda_[k] = params.row(k);}
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      { stat_lambda_[k].release();}
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      { stat_lambda_[k].resize(range);}
    }
};

/** @ingroup Clustering
 *  Specialization of the ParametersHandler struct for MixturePoisson_lk model
 **/
template <>
struct ParametersHandler<Clust::Poisson_lk_>
{
    typedef ModelParameters<Clust::Poisson_lk_> Parameters;
    /** Array of the lambdak_ statistics */
    Array1D< Stat::Online<Real, Real> > stat_lambda_;
    /** default constructor. All lambdas are initialized to 1. */
    inline ParametersHandler( int nbCluster)
                            : stat_lambda_(nbCluster)
    {}
    /** copy constructor.
     * @param param the parameters to copy.
     **/
    inline ParametersHandler( ParametersHandler const& param)
                            : stat_lambda_(param.stat_lambda_)
    {}
    /** destructor */
    inline ~ParametersHandler() {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    {
      stat_lambda_ = other.stat_lambda_;
      return *this;
    }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      { stat_lambda_[k].update(param.lambda_[k]);}
    }
    /** Set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      {
        param.lambda_[k] = stat_lambda_[k].mean();
        stat_lambda_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the lambda
     *  parameters on consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.lambda_.begin(); k<param.lambda_.end(); ++k)
      { params.row(k) = param.lambda_;}
    }
    /** Set the parameters of the mixture model.
     *  It is assumed that the array params store for each class the shapes and
     *  scales parameters on two consecutive rows.
     *  The number of column of params is the number of variables.
     **/
    template<class Array>
    inline void setParameters(Parameters& param, ExprBase<Array> const& params)
    {
      for(int k=param.lambda_.begin(); k<param.lambda_.end(); ++k)
      { param.lambda_[k] = params.row(k).mean();}
    }
    /** Release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      { stat_lambda_[k].release();}
    }
    /** Initialize the statistics. */
    inline void resize(Range const& range)
    {
      for(int k=stat_lambda_.begin(); k<stat_lambda_.end(); ++k)
      { stat_lambda_[k].release();}
    }
};

/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ Poisson mixture with the composer.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 */
template<int Id, class Data>
class PoissonBridge: public IMixtureBridge< PoissonBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< PoissonBridge<Id,Data> > Base;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits< PoissonBridge<Id,Data> >::ParamHandler ParamHandler;
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Poisson_
    };
    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::paramHandler_;
    using Base::p_data_;

    /** default constructor. Remove the missing values from the data set and
     *  initialize the mixture by setting the data set.
     *  @param p_data pointer on the DataBridge that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    PoissonBridge( DataBridge<Data>* p_data, std::string const& idData, int nbCluster)
                 : Base(p_data, idData, nbCluster)
    { removeMissing(); initializeBridge();}
    /** copy constructor */
    PoissonBridge( PoissonBridge const& bridge): Base(bridge)
    { initializeBridge();}
    /** destructor */
    virtual ~PoissonBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual PoissonBridge* clone() const { return new PoissonBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual PoissonBridge* create() const
    {
      PoissonBridge* p_bridge = new PoissonBridge( mixture_, this->idData(), this->nbCluster());
      p_bridge->p_data_ = p_data_;
      p_bridge->initializeBridge();
      return p_bridge;
    }
    /** This function is used in order to get the current values of the
     *  parameters.
     *  @param params the array with the parameters of the mixture.
     */
    template<class Array>
    void getParameters(Array& params) const;
    /** This function can be used to write summary of parameters to the output stream.
     *  @param out Stream where you want to write the summary of parameters.
     */
    virtual void writeParameters(std::ostream& out) const;

  private:
    /** This function will be used for the imputation of the missing data
     *  at the initialization.
     **/
    void removeMissing();
    /** This function will be used in order to initialize the mixture model
     *  using informations stored by the DataBridge. For example the missing
     *  values in the case of a DataBridge instance.
     **/
    void initializeBridge()
    {
      mixture_.setData(p_data_->dataij());
      paramHandler_.resize(p_data_->cols());
    }
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    PoissonBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                 : Base(mixture, idData, nbCluster)
    {}
    /** @return a safe value for the jth variable
     *  @param  dataij the matrix of the data set
     *  @param j index of the column with the safe value needed */
    static Type safeValue(Data const& dataij, int j)
    {

      int lmin = dataij.col(j).safe().minElt(), lmax = dataij.col(j).safe().maxElt();
      if (lmax -lmin > 10)
      { return Real(dataij.col(j).safe().sum())/dataij.sizeRows();}
      Array2DVector<int> count(Range(lmin, lmax, 0), 0);
      for (int i= dataij.beginRows(); i < dataij.endRows(); ++i)
      {
        if (!Arithmetic<int>::isNA(dataij(i,j)))
          count[dataij(i,j)]++;
      }
      int l; count.maxElt(l);
      return l;
    }

};

// implementation
template<int Id, class Data>
void PoissonBridge<Id, Data>::removeMissing()
{
  Type value = Type();
  int j, old_j = Arithmetic<int>::NA();
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  {
    j = it->second; // get column
    if (j != old_j)
    {
      old_j = j;
      value =  safeValue(p_data_->dataij(), j);
    }
    p_data_->dataij()(it->first, j) = value;
  }
}

template<int Id, class Data>
template<class Array>
void PoissonBridge<Id, Data>::getParameters(Array& params) const
{
  params.resize(this->nbCluster(), mixture_.p_data()->cols());
  for (int k= params.beginRows(); k < params.endRows(); ++k)
  {
    for (int j= mixture_.p_data()->beginCols();  j < mixture_.p_data()->endCols(); ++j)
    { params(k, j) = mixture_.lambda(k,j);}
  }
}

template<int Id, class Data>
void PoissonBridge<Id, Data>::writeParameters(ostream& os) const
{
  ArrayXX params;
  getParameters(params);
  for (int k= params.beginRows(); k < params.endRows(); ++k)
  {
    os << _T("---> Component ") << k << _T("\n");
    os << _T("lambda = ") << params.row(k);
  }
}

} // namespace STK

#endif /* STK_POISSONBRIDGE_H */
