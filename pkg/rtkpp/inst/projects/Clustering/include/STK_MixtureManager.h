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

namespace STK
{
/** @ingroup Clustering
 *  @brief A mixture manager is a factory class for injection dependency in the
 *  stk++ derived class of the IMixture interface
 *  (aka: the MixtureBridge<id> classes).
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
    // All Data type
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data Data_ajk_bjk;
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data Data_ajk_bj;
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data Data_sjk;
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data Data_sj;
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data Data_sk;
    typedef typename hidden::DataHandlerTraits<DataHandler, Real>::Data Data_s;
    typedef typename hidden::DataHandlerTraits<DataHandler, int>::Data Data_pjk;
    typedef typename hidden::DataHandlerTraits<DataHandler, int>::Data Data_pk;
    // All Data manager
    typedef DataManager<Clust::Gamma_ajk_bjk_, Data_ajk_bjk> DataManager_ajk_bjk;
    typedef DataManager<Clust::Gamma_ajk_bj_, Data_ajk_bj> DataManager_ajk_bj;
    typedef DataManager<Clust::Gaussian_sjk_, Data_sjk> DataManager_sjk;
    typedef DataManager<Clust::Gaussian_sj_, Data_sj> DataManager_sj;
    typedef DataManager<Clust::Gaussian_sk_, Data_sk> DataManager_sk;
    typedef DataManager<Clust::Gaussian_s_, Data_s> DataManager_s;
    typedef DataManager<Clust::Categorical_pjk_, Data_pjk> DataManager_pjk;
    typedef DataManager<Clust::Categorical_pk_, Data_pk> DataManager_pk;
    // All mixture bridges
    typedef MixtureBridge<Clust::Gamma_ajk_bjk_, Data_ajk_bjk> MixtureBridge_ajk_bjk;
    typedef MixtureBridge<Clust::Gamma_ajk_bj_, Data_ajk_bj> MixtureBridge_ajk_bj;
    typedef MixtureBridge<Clust::Gaussian_sjk_, Data_sjk> MixtureBridge_sjk;
    typedef MixtureBridge<Clust::Gaussian_sj_, Data_sj> MixtureBridge_sj;
    typedef MixtureBridge<Clust::Gaussian_sk_, Data_sk> MixtureBridge_sk;
    typedef MixtureBridge<Clust::Gaussian_s_, Data_s> MixtureBridge_s;
    typedef MixtureBridge<Clust::Categorical_pjk_, Data_pjk> MixtureBridge_pjk;
    typedef MixtureBridge<Clust::Categorical_pk_, Data_pk> MixtureBridge_pk;

    /** Default constructor, need an instance of a DataHandler.  */
    MixtureManager(DataHandler const& handler) : handler_(handler) {}
    /** destructor */
    ~MixtureManager()
    {
      for (DataIterator it = v_data_.begin() ; it != v_data_.end(); ++it)
      { delete (*it);}
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
     *  @param idModel id of the model
     *  @param idData name of the model
     *  @param nbCluster number of cluster of the model
     **/
    IMixture* createMixture(Clust::Mixture idModel, String const& idData, int nbCluster)
    {
      switch (idModel)
      {
        // gamma_ajk_bjk_ model
        case Clust::Gamma_ajk_bjk_:
        {
          DataManager_ajk_bjk* p_data = new DataManager_ajk_bjk();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_ajk_bjk* p_bridge = new MixtureBridge_ajk_bjk( p_data, idData, nbCluster);
          return p_bridge;
        }
        break;
        // gamma_ajk_bj_ model
        case Clust::Gamma_ajk_bj_:
        {
          DataManager_ajk_bj* p_data = new DataManager_ajk_bj();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_ajk_bj* p_bridge = new MixtureBridge_ajk_bj(p_data, idData, nbCluster);
          return p_bridge;
        }
        break;
        // Gaussian_sjk_ model
        case Clust::Gaussian_sjk_:
        {
          DataManager_sjk* p_data = new DataManager_sjk();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_sjk* p_bridge = new MixtureBridge_sjk(p_data, idData, nbCluster);
          stk_cout << p_bridge->idName() <<"\n";
          return p_bridge;
        }
        break;
        // Gaussian_sk_ model
        case Clust::Gaussian_sk_:
        {
          DataManager_sk* p_data = new DataManager_sk();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_sk* p_bridge = new MixtureBridge_sk(p_data, idData, nbCluster);
          return p_bridge;
        }
        break;
        // Gaussian_sj_ model
        case Clust::Gaussian_sj_:
        {
          DataManager_sj* p_data = new DataManager_sj();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_sj* p_bridge = new MixtureBridge_sj(p_data, idData, nbCluster);
          return p_bridge;
        }
        break;
        // Gaussian_s_ model
        case Clust::Gaussian_s_:
        {
          DataManager_s* p_data = new DataManager_s();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_s* p_bridge = new MixtureBridge_s(p_data, idData, nbCluster);
          return p_bridge;
        }
        break;
        // Categorical_pjk_ model
        case Clust::Categorical_pjk_:
        {
          DataManager_pjk* p_data = new DataManager_pjk();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_pjk* p_bridge = new MixtureBridge_pjk(p_data, idData, nbCluster);
          return p_bridge;
        }
        break;
        // Categorical_pjk_ model
        case Clust::Categorical_pk_:
        {
          DataManager_pk* p_data = new DataManager_pk();
          registerDataManager(p_data);
          handler_.getData(idData, p_data->m_dataij_, p_data->nbVariable_ );
          p_data->initialize();
          MixtureBridge_pk* p_bridge = new MixtureBridge_pk(p_data, idData, nbCluster);
          return p_bridge;
        }
        break;
        default:
          return 0; // 0 if idModel is not implemented
          break;
      }
      return 0; // 0 if idModel is not a stk++ model
    }
    /** get the parameters from an IMixture given its Id.*/
    void getParameters(IMixture* p_mixture, std::string idData, Array2D<Real>& data) const
    {
#ifdef STK_MIXTURE_VERY_VERBOSE
      stk_cout << _T("-------------------------------\n")
               << _T("Entering MixtureManager::getParameters with")
               << _T("idData = ") << idData << _T("\n");
#endif
      std::string idModelName;
      if (!handler_.getIdModel( idData, idModelName))
      {
#ifdef STK_MIXTURE_VERY_VERBOSE
        stk_cout << _T("In MixtureManager::getParameters, fail to get idData = ") << idData << _T("\n");
#endif
        return;
      }
#ifdef STK_MIXTURE_VERY_VERBOSE
      stk_cout << _T("In MixtureManager::getParameters, success to get idData = ") << idData << _T("\n");
      stk_cout << _T("In MixtureManager::getParameters, idModelName = ") << idModelName << _T("\n");
#endif
      Clust::Mixture idModel = Clust::stringToMixture(idModelName);
      // up-cast... (Yes it's bad....;)...)
      switch (idModel)
      {
        // gamma_ajk_bjk_ model
        case Clust::Gamma_ajk_bjk_:
        { static_cast<MixtureBridge_ajk_bjk const*>(p_mixture)->getParameters(data);}
        break;
        // gamma_ajk_bj_ model
        case Clust::Gamma_ajk_bj_:
        { static_cast<MixtureBridge_ajk_bj const*>(p_mixture)->getParameters(data);}
        break;
        // Gaussian_sjk_ model
        case Clust::Gaussian_sjk_:
        { static_cast<MixtureBridge_sjk const*>(p_mixture)->getParameters(data);}
        break;
        // Gaussian_sk_ model
        case Clust::Gaussian_sk_:
        { static_cast<MixtureBridge_sk const*>(p_mixture)->getParameters(data);}
        break;
        // Gaussian_sj_ model
        case Clust::Gaussian_sj_:
        { static_cast<MixtureBridge_sj const*>(p_mixture)->getParameters(data);}
        break;
        // Gaussian_s_ model
        case Clust::Gaussian_s_:
        { static_cast<MixtureBridge_s const*>(p_mixture)->getParameters(data);}
        break;
        // Categorical_pjk_ model
        case Clust::Categorical_pjk_:
        { static_cast<MixtureBridge_pjk const*>(p_mixture)->getParameters(data);}
        break;
        // Categorical_pjk_ model
        case Clust::Categorical_pk_:
        { static_cast<MixtureBridge_pk const*>(p_mixture)->getParameters(data);}
        break;
        default: // idModel is not implemented
        break;
      }
    }
  protected:
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
};

} // namespace STK

#endif /* STK_MIXTUREMANAGER_H */
