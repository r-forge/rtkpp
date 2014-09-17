/*--------------------------------------------------------------------*/
/*     Copyright (C) 2004-2014  Serge Iovleff

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

/** @file STK_MixtureManager.h
 *  @brief In this file we define the MixtureManager class.
 **/


#ifndef STK_MIXTUREMANAGER_H
#define STK_MIXTUREMANAGER_H

#include "MixturesBridges/STK_MixtureBridge.h"
#include "DManager/include/STK_IDataHandler.h"

#define STK_CREATE_MIXTURE(Data, Bridge) \
          Data* p_data = new Data(idData); \
          registerDataManager(p_data); \
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ ); \
          p_data->initialize(); \
          Bridge* p_bridge = new Bridge( p_data, idData, nbCluster);  \
          return p_bridge;

namespace STK
{
/** @ingroup Clustering
 *  @brief A mixture manager is a factory class for injection dependency in the
 *  stk++ derived class of the IMixture interface
 *  (aka: the MixtureBridge<Id> classes).
 *
 *  It handles all the creation and initialization stuff needed by the
 *  mixture models of the stkpp library.
 */
template<class DataHandler>
class MixtureManager
{
  public:
    typedef std::vector<IDataManager*>::const_iterator ConstDataIterator;
    typedef std::vector<IDataManager*>::iterator DataIterator;

    /** InfoMap is a map of pairs (idName, idModel) */
    typedef typename DataHandler::InfoMap InfoMap;
    // Real and Integer Data type
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data DataReal;
    typedef typename hidden::DataHandlerTraits<DataHandler, Integer>::Data DataInt;
    // All Data manager
    typedef DataManager<Clust::Gamma_, DataReal> DataManagerGamma;
    typedef DataManager<Clust::Gaussian_, DataReal> DataManagerGaussian;
    typedef DataManager<Clust::Categorical_, DataInt> DataManagerCategorical;
    // All Gamma bridges
    typedef MixtureBridge<Clust::Gamma_ajk_bjk_, DataReal> MixtureBridge_ajk_bjk;
    typedef MixtureBridge<Clust::Gamma_ajk_bk_,  DataReal> MixtureBridge_ajk_bk;
    typedef MixtureBridge<Clust::Gamma_ajk_bj_,  DataReal> MixtureBridge_ajk_bj;
    typedef MixtureBridge<Clust::Gamma_ajk_b_,   DataReal> MixtureBridge_ajk_b;
    typedef MixtureBridge<Clust::Gamma_ak_bjk_,  DataReal> MixtureBridge_ak_bjk;
    typedef MixtureBridge<Clust::Gamma_ak_bk_,   DataReal> MixtureBridge_ak_bk;
    typedef MixtureBridge<Clust::Gamma_ak_bj_,   DataReal> MixtureBridge_ak_bj;
    typedef MixtureBridge<Clust::Gamma_ak_b_,    DataReal> MixtureBridge_ak_b;
    typedef MixtureBridge<Clust::Gamma_aj_bjk_,  DataReal> MixtureBridge_aj_bjk;
    typedef MixtureBridge<Clust::Gamma_aj_bk_,   DataReal> MixtureBridge_aj_bk;
    typedef MixtureBridge<Clust::Gamma_a_bjk_,   DataReal> MixtureBridge_a_bjk;
    typedef MixtureBridge<Clust::Gamma_a_bk_,    DataReal> MixtureBridge_a_bk;
    // All Gaussian bridges
    typedef MixtureBridge<Clust::Gaussian_sjk_, DataReal> MixtureBridge_sjk;
    typedef MixtureBridge<Clust::Gaussian_sj_,  DataReal> MixtureBridge_sj;
    typedef MixtureBridge<Clust::Gaussian_sk_,  DataReal> MixtureBridge_sk;
    typedef MixtureBridge<Clust::Gaussian_s_,   DataReal> MixtureBridge_s;
    // All Categorical bridges
    typedef MixtureBridge<Clust::Categorical_pjk_, DataInt> MixtureBridge_pjk;
    typedef MixtureBridge<Clust::Categorical_pk_,  DataInt> MixtureBridge_pk;

    /** Default constructor, need an instance of a DataHandler.  */
    MixtureManager(DataHandler const& handler) : handler_(handler) {}
    /** destructor */
    ~MixtureManager()
    {
      for (DataIterator it = v_data_.begin() ; it != v_data_.end(); ++it)
      { delete (*it);}
    }
    /** Utility function allowing to find the idModel from the idData
     *  @param idData the id name of the data we want the idModel
     *  @return the idModel
     **/
    Clust::Mixture getIdModel( String const& idData) const
    {
      std::string idModelName;
      if (!handler_.getIdModel( idData, idModelName))
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("In MixtureManager::getIdModel, fail to get idData = ") << idData << _T("\n");
#endif
        return Clust::unknown_mixture_;
      }
#ifdef STK_MIXTURE_VERY_VERBOSE
      stk_cout << _T("In MixtureManager::getIdModel, success to get idData = ") << idData << _T("\n");
      stk_cout << _T("In MixtureManager::getIdModel, idModelName = ") << idModelName << _T("\n");
#endif
      return Clust::stringToMixture(idModelName);
    }
    /** Utility function allowing to create and register all the stk++ mixtures
     *  defined in the handler.
     *  @param composer the composer claiming the mixtures
     *  @param nbCluster the number of clusters
     **/
    void createMixtures(MixtureComposer& composer, int nbCluster)
    {
      for (typename InfoMap::const_iterator it=handler_.info().begin(); it!=handler_.info().end(); ++it)
      {
        std::string idData = it->first;
        Clust::Mixture idModel = Clust::stringToMixture(it->second);
        // get a mixture
        IMixture* p_bridge = createMixture(idModel, idData, nbCluster);
        if (p_bridge) composer.registerMixture(p_bridge);
      }
    }
    /** create a mixture and initialize it.
     *  @param idData name of the model
     *  @param nbCluster number of cluster of the model
     *  @return 0 if the idData is not find, the result of @c createMixture( idModel, idData, nbCluster)
     *  otherwise.
     **/
    IMixture* createMixture(String const& idData, int nbCluster)
    {
      std::string idModelName;
      if (!handler_.getIdModel( idData, idModelName)) { return 0;};
      Clust::Mixture idModel = Clust::stringToMixture(idModelName);
      return createMixture( idModel, idData, nbCluster);
    }
    /** create a mixture and initialize it.
     *  @param idModel Id name of the model
     *  @param idData Id name of the data
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixture(Clust::Mixture idModel, String const& idData, int nbCluster)
    {
      switch (idModel)
      {
        // gamma_ajk_bjk_ model
        case Clust::Gamma_ajk_bjk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ajk_bjk)}
        break;
        // gamma_ajk_bk_ model
        case Clust::Gamma_ajk_bk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ajk_bk)}
        break;
        // gamma_ajk_bj_ model
        case Clust::Gamma_ajk_bj_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ajk_bj)}
        break;
        // gamma_ajk_b_ model
        case Clust::Gamma_ajk_b_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ajk_b)}
        break;
        // gamma_ak_bjk_ model
        case Clust::Gamma_ak_bjk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ak_bjk)}
        break;
        // gamma_ak_bk_ model
        case Clust::Gamma_ak_bk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ak_bk)}
        break;
        // gamma_ak_bj_ model
        case Clust::Gamma_ak_bj_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ak_bj)}
        break;
        // gamma_ajk_b_ model
        case Clust::Gamma_ak_b_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_ak_b)}
        break;
        // gamma_aj_bjk_ model
        case Clust::Gamma_aj_bjk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_aj_bjk)}
        break;
        // gamma_aj_bk_ model
        case Clust::Gamma_aj_bk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_aj_bk)}
        // gamma_aj_bjk_ model
        case Clust::Gamma_a_bjk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_a_bjk)}
        break;
        // gamma_aj_bk_ model
        case Clust::Gamma_a_bk_:
        { STK_CREATE_MIXTURE(DataManagerGamma, MixtureBridge_a_bk)}
        // Gaussian_sjk_ model
        case Clust::Gaussian_sjk_:
        { STK_CREATE_MIXTURE(DataManagerGaussian, MixtureBridge_sjk)}
        break;
        // Gaussian_sk_ model
        case Clust::Gaussian_sk_:
        { STK_CREATE_MIXTURE(DataManagerGaussian, MixtureBridge_sk)}
        break;
        // Gaussian_sj_ model
        case Clust::Gaussian_sj_:
        { STK_CREATE_MIXTURE(DataManagerGaussian, MixtureBridge_sj)}
        break;
        // Gaussian_s_ model
        case Clust::Gaussian_s_:
        { STK_CREATE_MIXTURE(DataManagerGaussian, MixtureBridge_s)}
        break;
        // Categorical_pjk_ model
        case Clust::Categorical_pjk_:
        { STK_CREATE_MIXTURE(DataManagerCategorical, MixtureBridge_pjk)}
        break;
        // Categorical_pjk_ model
        case Clust::Categorical_pk_:
        { STK_CREATE_MIXTURE(DataManagerCategorical, MixtureBridge_pk)}
        break;
        default:
          return 0; // 0 if idModel is not implemented
          break;
      }
      return 0; // 0 if idModel is not a stk++ model
    }
    /** release a data set from v_data_.
     *  @param idData name of the data set to release
     **/
    void releaseMixture(String const& idData)
    {
      for (DataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
      { if ((*it)->idData() == idData) {delete (*it); v_data_.erase(it); break;}}
    }
    /** get the parameters from an IMixture.
     *  @param p_mixture pointer on the mixture
     *  @param idData Id name of the data set attached to the mixture
     *  @param data the array to return with the parameters
     **/
    void getParameters(IMixture* p_mixture, std::string idData, Array2D<Real>& data) const
    {
      Clust::Mixture idModel = getIdModel(idData);
      if (idModel == Clust::unknown_mixture_) return;
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // gamma models
        case Clust::Gamma_ajk_bjk_:
        { static_cast<MixtureBridge_ajk_bjk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_ajk_bk_:
        { static_cast<MixtureBridge_ajk_bk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_ajk_bj_:
        { static_cast<MixtureBridge_ajk_bj const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_ajk_b_:
        { static_cast<MixtureBridge_ajk_b const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_ak_bjk_:
        { static_cast<MixtureBridge_ak_bjk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_ak_bk_:
        { static_cast<MixtureBridge_ak_bk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_ak_bj_:
        { static_cast<MixtureBridge_ak_bj const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_ak_b_:
        { static_cast<MixtureBridge_ak_b const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_aj_bjk_:
        { static_cast<MixtureBridge_aj_bjk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_aj_bk_:
        { static_cast<MixtureBridge_aj_bk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_a_bjk_:
        { static_cast<MixtureBridge_a_bjk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gamma_a_bk_:
        { static_cast<MixtureBridge_a_bk const*>(p_mixture)->getParameters(data);}
        break;
        // Gaussian models
        case Clust::Gaussian_sjk_:
        { static_cast<MixtureBridge_sjk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gaussian_sk_:
        { static_cast<MixtureBridge_sk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gaussian_sj_:
        { static_cast<MixtureBridge_sj const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Gaussian_s_:
        { static_cast<MixtureBridge_s const*>(p_mixture)->getParameters(data);}
        break;
        // Categorical models
        case Clust::Categorical_pjk_:
        { static_cast<MixtureBridge_pjk const*>(p_mixture)->getParameters(data);}
        break;
        case Clust::Categorical_pk_:
        { static_cast<MixtureBridge_pk const*>(p_mixture)->getParameters(data);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
    /** get the missing values of a data set.
     *  @param idData Id name of the data set attached to the mixture
     *  @param data the array to return with the missing values
     **/
    template<typename Type>
    void getMissingValues( String const& idData, std::vector< std::pair< std::pair<int,int>, Type > >& data) const
    {
      Clust::Mixture idModel = getIdModel(idData);
      if (idModel == Clust::unknown_mixture_) return;
      Clust::MixtureClass idClass = Clust::MixtureToMixtureClass(idModel);
      IDataManager* p_manager = getDataManager(idData);
      // up-cast... (Yes it's bad....;)...)
      switch (idClass)
      {
        // gamma models
        case Clust::Gamma_:
        { static_cast<DataManagerGamma const*>(p_manager)->getMissingValues(data);}
        break;
        // Gaussian_ models
        case Clust::Gaussian_:
        { static_cast<DataManagerGaussian const*>(p_manager)->getMissingValues(data);}
        break;
        // Categorical_ models
        case Clust::Categorical_:
        { static_cast<DataManagerCategorical const*>(p_manager)->getMissingValues(data);}
        break;
        default: // idClass is not implemented
        break;
      }
    }

    /** @brief register a data manager to the MixtureManager.
     *  For each mixture created and registered, a data manager is created
     *  and registered so that it will be deleted when the mixture itself is
     *  deleted.
     *  @param p_data a pointer on the data manager
     **/
    inline void registerDataManager(IDataManager* p_data)
    { v_data_.push_back(p_data);}

  private:
    /** pointer to the dataHandler */
    DataHandler const& handler_;
    /** vector of pointers to the data components */
    std::vector<IDataManager*> v_data_;
    /** Utility lookup function allowing to find a DataManager from its idData
     *  @param idData the id name of the mixture we want to get
     *  @return a pointer on the DataManager
     **/
    IDataManager* getDataManager( String const& idData) const
    {
      for (ConstDataIterator it = v_data_.begin(); it != v_data_.end(); ++it)
      {  if ((*it)->idData() == idData) return (*it);}
      return 0;
    }
};

} // namespace STK

#undef STK_CREATE_MIXTURE

#endif /* STK_MIXTUREMANAGER_H */
