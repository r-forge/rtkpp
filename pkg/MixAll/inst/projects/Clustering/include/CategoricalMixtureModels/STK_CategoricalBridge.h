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

/** @file STK_CategoricalBridge.h
 *  @brief In this file we define the bridge classes between the categorical
 *  mixtures and the composer.
 **/

#ifndef STK_CATEGORICALBRIDGE_H
#define STK_CATEGORICALBRIDGE_H

#include "STK_MixtureCategorical_pk.h"
#include "STK_MixtureCategorical_pjk.h"

#include "../STK_IMixtureBridge.h"

#include <DManager/include/STK_DataBridge.h>
#include <STatistiK/include/STK_Stat_Online.h>

namespace STK
{

// forward declaration
template<int Id, class Data> class CategoricalBridge;

namespace hidden
{
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixtureCategorical_pjk model
 **/
template<class Data_>
struct MixtureBridgeTraits< CategoricalBridge<Clust::Categorical_pjk_, Data_> >
{
  /** Structure storing the Data */
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the Mixture model */
  typedef MixtureCategorical_pjk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Categorical_pjk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Categorical_pjk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Categorical_
  };
};
/** @ingroup hidden
 *  Partial specialization of the MixtureBridgeTraits for the MixtureCategorical_pk model
 **/
template<class Data_>
struct MixtureBridgeTraits< CategoricalBridge< Clust::Categorical_pk_, Data_> >
{
  /** Structure storing the Data */
  typedef Data_ Data;
  /** Data Type */
  typedef typename Data_::Type Type;
  /** Type of the mixture model */
  typedef MixtureCategorical_pk<Data> Mixture;
  /** Type of the structure storing the mixture parameters */
  typedef ModelParameters<Clust::Categorical_pk_> Parameters;
  /** Type of the parameter handler */
  typedef ParametersHandler<Clust::Categorical_pk_> ParamHandler;
  enum
  {
    idMixtureClass_ = Clust::Categorical_
  };
};

} // namespace hidden

/** @ingroup Clustering
 * Specialization of the ParametersHandler struct for MixtureCategorical_pjk model */
template <>
struct ParametersHandler<Clust::Categorical_pjk_>
{
    typedef ModelParameters<Clust::Categorical_pjk_> Parameters;
    /** statistics of the probabilities */
    Array1D<  Stat::Online<CArrayXX, Real>  > stat_proba_;
    /** default constructor */
    ParametersHandler(int nbCluster): stat_proba_(nbCluster) {}
    /** copy constructor */
    ParametersHandler(ParametersHandler const& model): stat_proba_(model.stat_proba_) {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    { stat_proba_ = other.stat_proba_; return *this; }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      { stat_proba_[k].update(param.proba_[k]);}
    }
    /** set and release the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      {
        param.proba_[k] = stat_proba_[k].mean();
        stat_proba_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the array of
     *  probabilities on consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.proba_.begin(), kp= params.beginRows(); k<param.proba_.end(); ++k)
      {
        for (int l = param.proba_[k].beginRows(); l < param.proba_[k].endRows(); ++l, ++kp)
        {
          for (int j = param.proba_[k].beginCols(); j < param.proba_[k].endCols(); ++j)
          { params(kp , j) = param.proba_[k](l, j);}
        }
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
      int kp = params.beginRows();
      for(int k=param.proba_.begin(); k<param.proba_.end(); ++k)
      {
        for (int l = param.proba_[k].beginRows(); l < param.proba_[k].endRows(); ++l, ++kp)
        {
          for (int j = param.proba_[k].beginCols(); j < param.proba_[k].endCols(); ++j)
          { param.proba_[k](l, j) = params(kp , j) ;}
        }
      }
    }
    /** Set the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      { stat_proba_[k].release();}
    }
    /** Initialize the parameters of the model. */
    void resize(Range const& rangeModalities, Range const& rangeData )
    {
      for (int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      {
        stat_proba_[k].resize(rangeModalities, rangeData);
      }
    }
};

/** Specialization of the ParametersHandler struct for MixtureCategorical_pk model */
template <>
struct ParametersHandler<Clust::Categorical_pk_>
{
    typedef ModelParameters<Clust::Categorical_pk_> Parameters;
    /** statistics of the probabilities */
    Array1D<  Stat::Online<CVectorX, Real>  > stat_proba_;
    /** default constructor */
    ParametersHandler(int nbCluster): stat_proba_(nbCluster) {}
    /** copy constructor */
    ParametersHandler(ParametersHandler const& model): stat_proba_(model.stat_proba_) {}
    /** copy operator */
    inline ParametersHandler& operator=( ParametersHandler const& other)
    { stat_proba_ = other.stat_proba_; return *this; }
    /** update statistics of the parameters. */
    inline void updateStatistics(Parameters const& param)
    {
      for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      { stat_proba_[k].update(param.proba_[k]);}
    }
    /** set the computed statistics */
    inline void setStatistics(Parameters& param)
    {
      for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      {
        param.proba_[k] = stat_proba_[k].mean();
        stat_proba_[k].release();
      }
    }
    /** Get the parameters of the mixture model.
     *  It is assumed that the array params store for each class the array of
     *  probabilities on consecutive rows.
     *  The number of column of params is the number of variables.
     *  @note the array params has to be resized before any call
     **/
    template<class Array>
    inline void getParameters(Parameters const& param, ArrayBase<Array>& params)
    {
      for(int k=param.proba_.begin(), kp= params.beginRows(); k<param.proba_.end(); ++k)
      {
        for (int l = param.proba_[k].beginRows(); l < param.proba_[k].endRows(); ++l, ++kp)
        {
          for (int j = param.proba_[k].beginCols(); j < param.proba_[k].endCols(); ++j)
          { params(kp , j) = param.proba_[k][l];}
        }
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
      int kp = params.beginRows();
      for(int k=param.proba_.begin(); k<param.proba_.end(); ++k)
      {
        for (int l = param.proba_[k].beginRows(); l < param.proba_[k].endRows(); ++l, ++kp)
        {
          param.proba_[k][l] = 0.;
          for (int j = params.beginCols(); j < params.endCols(); ++j)
          { param.proba_[k][l] += params(kp , j) ;}
          param.proba_[k][l] /= param.proba_[k].sizeCols();
        }
      }
    }
    /** release the computed statistics */
    inline void releaseStatistics()
    {
      for(int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      { stat_proba_[k].release();}
    }
    /** Initialize the parameters of the model. */
    void resize(Range const& rangeModalities, Range const& rangeData )
    {
      for (int k=stat_proba_.begin(); k<stat_proba_.end(); ++k)
      { stat_proba_[k].resize(rangeModalities);}
    }
};

/** @ingroup Clustering
 *  @brief Templated implementation of the IMixture interface allowing
 *  to bridge a STK++ mixture with the composer.
 *
 *  This class inherit from the IMixtureBridge.
 *
 * @tparam Id is any identifier of a concrete model deriving from the
 * interface STK::IMixtureDensity class.
 */
template<int Id, class Data>
class CategoricalBridge: public IMixtureBridge< CategoricalBridge<Id,Data> >
{
  public:
    // Base class
    typedef IMixtureBridge< CategoricalBridge<Id,Data> > Base;
    // type of Mixture
    typedef typename hidden::MixtureBridgeTraits< CategoricalBridge<Id,Data> >::Mixture Mixture;
    typedef typename hidden::MixtureBridgeTraits< CategoricalBridge<Id,Data> >::Parameters Parameters;
    typedef typename hidden::MixtureBridgeTraits< CategoricalBridge<Id,Data> >::ParamHandler ParamHandler;
    // type of data
    typedef typename Data::Type Type;
    // class of mixture
    enum
    {
      idMixtureClass_ = Clust::Categorical_
    };

    typedef std::vector<std::pair<int,int> >::const_iterator ConstIterator;
    using Base::mixture_;
    using Base::paramHandler_;
    using Base::p_data_;
    using Base::p_tik;

    /** default constructor.
     *  - Remove the missing values from the data set,
     *  - initialize the mixture by setting the data set
     *  - initialize the ParamHandler
     *  @param p_data pointer on the DataBridge that will be used by the bridge.
     *  @param idData id name of the mixture model
     *  @param nbCluster number of cluster
     **/
    CategoricalBridge( DataBridge<Data>* p_data, std::string const& idData, int nbCluster)
                     : Base( p_data, idData, nbCluster)
    { removeMissing(); initializeBridge();}
    /** copy constructor */
    CategoricalBridge( CategoricalBridge const& bridge): Base(bridge)
    { initializeBridge();}
    /** destructor */
    virtual ~CategoricalBridge() {}
    /** This is a standard clone function in usual sense. It must be defined to
     *  provide new object of your class with values of various parameters
     *  equal to the values of calling object. In other words, this is
     *  equivalent to polymorphic copy constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual CategoricalBridge* clone() const { return new CategoricalBridge(*this);}
    /** This is a standard create function in usual sense. It must be defined to
     *  provide new object of your class with correct dimensions and state.
     *  In other words, this is equivalent to virtual constructor.
     *  @return New instance of class as that of calling object.
     */
    virtual CategoricalBridge* create() const
    {
      CategoricalBridge* p_bridge = new CategoricalBridge( mixture_, this->idData(), this->nbCluster());
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
    /** Write the parameters on the output stream os */
    void writeParameters(ostream& os) const;

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
      paramHandler_.resize(mixture_.modalities(), p_data_->cols());
    }
    /** protected constructor to use in order to create a bridge.
     *  @param mixture the mixture to copy
     *  @param idData id name of the mixture
     *  @param nbCluster number of cluster
     **/
    CategoricalBridge( Mixture const& mixture, std::string const& idData, int nbCluster)
                     : Base(mixture, idData, nbCluster)
    {}
    /** @return a safe value for the jth variable
     *  @param  dataij the matrix of the data set
     *  @param j index of the column with the safe value needed */
    Type safeValue(Data const& dataij, int j)
    {
       int lmin = dataij.col(j).safe().minElt(), lmax = dataij.col(j).safe().maxElt();
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
void CategoricalBridge<Id, Data>::removeMissing()
{
  Type value = Type();
  int j, old_j = Arithmetic<int>::NA();
  for(ConstIterator it = p_data_->v_missing().begin(); it!= p_data_->v_missing().end(); ++it)
  {
    j = it->second; // get column
    if (j != old_j)
    {
      old_j = j;
      value =  safeValue(p_data_->dataij(),j);
    }
    p_data_->dataij()(it->first, j) = value;
  }
}

template<int Id, class Data>
template<class Array>
void CategoricalBridge<Id, Data>::getParameters(Array& params) const
{
    int nbCluster    = this->nbCluster();
    int nbModalities = mixture_.modalities().size();

    params.resize(nbModalities * nbCluster, mixture_.p_data()->cols());
    for (int k = 0; k < nbCluster; ++k)
    {
      for (int j = mixture_.p_data()->beginCols(); j < mixture_.p_data()->endCols(); ++j)
      {
        for (int l = 0; l < nbModalities; ++l)
        { params(baseIdx+k * nbModalities + l, j) = mixture_.proba(baseIdx+k, j, mixture_.modalities().begin() + l);}
      }
    }
}

/* Write the parameters on the output stream os */
template<int Id, class Data>
void CategoricalBridge<Id, Data>::writeParameters(ostream& os) const
{
  ArrayXX p(mixture_.modalities(), mixture_.p_data()->cols());
  for (int k= p_tik()->beginCols(); k < p_tik()->endCols(); ++k)
  {
    // store proba values in an array for a nice output
    for (int j= p.beginCols();  j < p.endCols(); ++j)
    {
      for (int l= mixture_.modalities().begin(); l < mixture_.modalities().end(); ++l)
      { p(l, j) = mixture_.proba(k,j,l);}
    }
    os << _T("---> Component ") << k << _T("\n");
    os << _T("probabilities =\n") << p  << _T("\n");
  }
}

} // namespace STK

#endif /* STK_CATEGORICALBRIDGE_H */
