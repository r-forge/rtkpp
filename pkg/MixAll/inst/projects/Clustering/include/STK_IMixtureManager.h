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
 * created on: 15 mars 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_IMixtureManager.h
 *  @brief In this file we define the Interface IMixtureManager class.
 **/


#ifndef STK_IMIXTUREMANAGER_H
#define STK_IMIXTUREMANAGER_H

#include <DManager/include/STK_DataHandlerBase.h>
#include <DManager/include/STK_DataBridge.h>

#include <Arrays/include/STK_Array2D.h> // for get and set parameters

#include "STK_IMixture.h"

namespace STK
{
/** @ingroup Clustering
 *  @brief Interface base class for mixture managers.
 *
 *  A mixture manager is a factory class for injection dependency in the
 *  STK++ derived class of the IMixtureComposer:
 *  - It handles all the creation and initialization stuff needed by mixture models,
 *  - It allows to get parameters and imputed missing values from specific mixtures,
 *  - It allows also to set parameters to a specific mixture,
 *  - all data set are enclosed in a DataBridge structure and stored in vector v_data_
 *
 *  The pure virtual method to implement in derived classes are
 *  @code
 *    virtual void getParameters(IMixture* p_mixture, ArrayXX& data) const =0;
 *    virtual void setParameters(IMixture* p_mixture, ArrayXX const& data) const =0;
 *    virtual IMixture* createMixtureImpl(String const& modelName, String const& idData, int nbCluster) =0;
 *  @endcode
 *
 *  @tparam DataHandler any concrete class from the interface STK::DataHandlerBase
 */
template<class DataHandler>
class IMixtureManager
{
  public:
    /** Default constructor, need an instance of a DataHandler. */
    IMixtureManager(DataHandler const* const p_handler);
    /** destructor */
    virtual ~IMixtureManager();
    /** @return constant pointer on the data handler */
    DataHandler const* const p_handler() const { return p_handler_;}

    /** Utility function allowing to find the idModel from the idData
     *  @param idData the id name of the data we want the idModel
     *  @return the idModel
     **/
    Clust::Mixture getIdModel( String const& idData) const;
    /** Utility function allowing to find the idModel name from the idData
     *  @param idData the id name of the data we want the idModel
     *  @return the idModel name
     **/
    String getIdModelName( String const& idData) const;
    /** @brief create a mixture and initialize it.
     *  This method get the modelName from the DataHandler and then delegate
     *  the concrete creation to derived class using the pure virtual method
     *   @c createMixtureImpl.
     *  @param idData name of the model
     *  @param nbCluster number of cluster of the model
     *  @return 0 if the idData is not find, the result of
     *  @c createMixture( modelName, idData, nbCluster) otherwise.
     **/
    IMixture* createMixture(String const& idData, int nbCluster);
    /** @brief register a data bridge to the IMixtureManager.
     *  For each mixture created and registered, a data manager is created
     *  and registered so that it will be deleted when the mixture itself is
     *  deleted.
     *  @param p_data a pointer on the data manager
     **/
    void registerDataBridge(IDataBridge* p_data);
    /** release a data bridge from v_data_.
     *  @param idData name of the data set to release
     **/
    void releaseDataBridge(String const& idData);

    // pure virtual methods
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param data the array to return with the parameters
     **/
    virtual void getParameters(IMixture* p_mixture, ArrayXX& data) const =0;
    /** set the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param data the array with the parameters to set
     **/
    virtual void setParameters(IMixture* p_mixture, ArrayXX const& data) const =0;

    /** get the wrapper for any kind of data set using its Id
     *  @param idData Id name of the data set attached to the mixture
     *  @return a constant reference on the array with the data set
     **/
    template<typename Type>
    typename hidden::DataHandlerTraits<DataHandler, Type>::Data const& getData( String const& idData) const;

  protected:
    /** Utility lookup function allowing to find a DataBridge from its idData
     *  @param idData the id name of the mixture we want to get
     *  @return a pointer on the DataBridge
     **/
    IDataBridge* getDataBridge( String const& idData) const;

  private:
    /** create a concrete mixture and initialize it.
     *  @param modelName, idData strings with the Id name of the model and of the data
     *  @param nbCluster number of cluster of the model
     **/
    virtual IMixture* createMixtureImpl(String const& modelName, String const& idData, int nbCluster) =0;
    /** A pointer on the concrete instance of the data handler */
    DataHandler const* const p_handler_;
    /** vector of pointers to the data components */
    std::vector<IDataBridge*> v_data_;
};

/* Utility function allowing to find the idModel from the idData
 *  @param idData the id name of the data we want the idModel
 *  @return the idModel
 **/
template<class DataHandler>
Clust::Mixture IMixtureManager<DataHandler>::getIdModel( String const& idData) const
{
  std::string modelName;
  if (!p_handler()->getIdModelName( idData, modelName))
  {
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("In IMixtureManager::getIdModel, fail to get idData = ") << idData << _T("\n");
#endif
    return Clust::unknown_mixture_;
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In IMixtureManager::getIdModel, success to get idData = ") << idData << _T("\n");
  stk_cout << _T("In IMixtureManager::getIdModel, modelName = ") << modelName << _T("\n");
#endif
  return Clust::stringToMixture(modelName);
}

/* Default constructor, need an instance of a DataHandler.  */
template<class DataHandler>
IMixtureManager<DataHandler>::IMixtureManager( DataHandler const* const p_handler)
                                             : p_handler_(p_handler)
{}
/* destructor */
template<class DataHandler>
IMixtureManager<DataHandler>::~IMixtureManager()
{
  typedef std::vector<IDataBridge*>::iterator DataIterator;
  for (DataIterator it = v_data_.begin() ; it != v_data_.end(); ++it)
  { delete (*it);}
}
/* Utility function allowing to find the idModel name from the idData
 *  @param idData the id name of the data we want the idModel
 *  @return the idModel name
 **/
template<class DataHandler>
String IMixtureManager<DataHandler>::getIdModelName( String const& idData) const
{
  std::string modelName;
  if (!p_handler_->getIdModelName( idData, modelName))
  {
#ifdef STK_MIXTURE_VERY_VERBOSE
    stk_cout << _T("In IMixtureManager::getIdModelName, fail to get idData = ") << idData << _T("\n");
#endif
  }
#ifdef STK_MIXTURE_VERY_VERBOSE
  stk_cout << _T("In IMixtureManager::getIdModeName, success to get idData = ") << idData << _T("\n");
  stk_cout << _T("In IMixtureManager::getIdModel, modelName = ") << modelName << _T("\n");
#endif
  return modelName;
}
/* create a mixture and initialize it.
 *  @param idData name of the model
 *  @param nbCluster number of cluster of the model
 *  @return 0 if the idData is not find, the result of @c createMixture( idModel, idData, nbCluster)
 *  otherwise.
 **/
template<class DataHandler>
IMixture* IMixtureManager<DataHandler>::createMixture(String const& idData, int nbCluster)
{
  std::string modelName;
  if (!p_handler_->getIdModelName( idData, modelName)) { return 0;};
  return createMixtureImpl( modelName, idData, nbCluster);
}
/* @brief register a data manager to the IMixtureManager.
 *  For each mixture created and registered, a data manager is created
 *  and registered so that it will be deleted when the mixture itself is
 *  deleted.
 *  @param p_data a pointer on the data manager
 **/
template<class DataHandler>
void IMixtureManager<DataHandler>::registerDataBridge(IDataBridge* p_data)
{ v_data_.push_back(p_data);}
/* release a data set from v_data_.
 *  @param idData name of the data set to release
 **/
template<class DataHandler>
void IMixtureManager<DataHandler>::releaseDataBridge(String const& idData)
{
  typedef std::vector<IDataBridge*>::iterator DataIterator;
  for (DataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  { if ((*it)->idData() == idData) {delete (*it); v_data_.erase(it); break;}}
}

/* Utility lookup function allowing to find a DataBridge from its idData
 *  @param idData the id name of the mixture we want to get
 *  @return a pointer on the DataBridge
 **/
template<class DataHandler>
IDataBridge* IMixtureManager<DataHandler>::getDataBridge( String const& idData) const
{
  typedef std::vector<IDataBridge*>::const_iterator ConstDataIterator;
  for (ConstDataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
  {  if ((*it)->idData() == idData) return (*it);}
  return 0;
}

/* get the wrapper for any kind of data set using its Id
 *  @param idData Id name of the data set attached to the mixture
 *  @return a constant reference on the array with the data set
 **/
template<class DataHandler>
template<typename Type>
typename hidden::DataHandlerTraits<DataHandler, Type>::Data const& IMixtureManager<DataHandler>::getData( String const& idData) const
{
  typedef typename hidden::DataHandlerTraits<DataHandler, Type>::Data DataType;
  typedef DataBridge<DataType> DataBridgeType;
  return static_cast<DataBridgeType*>(getDataBridge(idData))->dataij();
}
} // namespace STK


#endif /* STK_IMIXTUREMANAGER_H */
