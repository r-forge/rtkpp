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
 **/

/** @file STK_DataBridge.h
 *  @brief In this file we define the data wrapper class
 **/

#ifndef STK_DATABRIDGE_H
#define STK_DATABRIDGE_H

#include "STK_IDataBridge.h"

#include <Arrays/include/STK_ITContainer2D.h>

namespace STK
{
// forward declaration
template<class Data> class DataBridge;
namespace hidden
{
/** @ingroup hidden
 *  @brief Specialization of the Traits class for CArray class.
 */
template<class Data>
struct Traits< DataBridge<Data> >
{
  public:
    typedef typename Traits<Data>::Number Number;
    typedef typename Traits<Data>::Row Row;
    typedef typename Traits<Data>::Col Col;
    typedef typename Traits<Data>::SubRow SubRow;
    typedef typename Traits<Data>::SubCol SubCol;
    typedef typename Traits<Data>::SubVector SubVector;
    typedef typename Traits<Data>::SubArray SubArray;

    typedef typename Traits<Data>::Type Type;
    typedef typename Traits<Data>::ConstReturnType ConstReturnType;
    enum
    {
      structure_ = Traits<Data>::structure_,
      orient_    = Traits<Data>::orient_,
      sizeRows_  = Traits<Data>::sizeRows_,
      sizeCols_  = Traits<Data>::sizeCols_,
      storage_   = Traits<Data>::storage_
    };
};
} // hidden

/** @ingroup Clustering
 *  @brief bridge a data set in order to handle its missing values.
 *
 * @tparam Data The data bridged by the DataBridge class
 */
template<class Data>
class DataBridge: public IDataBridge, public ITContainer2D< DataBridge<Data> >
{
  public:
    typedef IDataBridge IBase;
    typedef ITContainer2D< DataBridge<Data> > Base;
    typedef std::vector<std::pair<int,int> > MissingIndexes;

    typedef typename hidden::Traits< DataBridge<Data> >::Type Type;
    typedef typename hidden::Traits< DataBridge<Data> >::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits< DataBridge<Data> >::Row Row;
    typedef typename hidden::Traits< DataBridge<Data> >::Col Col;
    typedef typename hidden::Traits< DataBridge<Data> >::SubRow SubRow;
    typedef typename hidden::Traits< DataBridge<Data> >::SubCol SubCol;
    typedef typename hidden::Traits< DataBridge<Data> >::SubArray SubArray;
    typedef typename hidden::Traits< DataBridge<Data> >::SubVector SubVector;

    enum
    {
      structure_ = hidden::Traits< DataBridge<Data> >::structure_,
      orient_    = hidden::Traits< DataBridge<Data> >::orient_,
      sizeRows_  = hidden::Traits< DataBridge<Data> >::sizeRows_,
      sizeCols_  = hidden::Traits< DataBridge<Data> >::sizeCols_,
      storage_   = hidden::Traits< DataBridge<Data> >::storage_
    };

    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    using IBase::v_missing_;

    /** default constructor. */
    inline DataBridge( std::string const& idData)
                     : IBase(idData), Base(), dataij_(), nbVariable_(0) {}
    /** constructor with data. */
    inline DataBridge( std::string const& idData, Data const& dataij)
                     : IBase(idData)
                     , Base(dataij.rows(), dataij.cols())
                     , dataij_(dataij)
                     , nbVariable_(dataij.sizeCols()) {}
    /** copy constructor (Warning: will copy the data set)
     *  @param bridge the DataBridge to copy
     **/
    DataBridge( DataBridge const& bridge)
              : IBase(bridge), Base(bridge)
              , dataij_(bridge.dataij_)
              , nbVariable_(bridge.nbVariable_) {}
    /** destructor */
    virtual ~DataBridge() {}
    /** getter. @return the number of variables (the number of columns of the data)  */
    inline int nbVariable() const { return nbVariable_;}
    /** getter. @return a constant reference on the data set*/
    inline Data const& dataij() const { return dataij_;}

    /** setter. @return a reference on the number of variables */
    inline int& nbVariable() { return nbVariable_;}
    /** setter. @return a reference on the data set
     *  @warning call to initialize is required if dataij_ is modified
     **/
    inline Data& dataij() { return dataij_;}

    /** cast operator */
    inline operator Data() const { return dataij_;}

   /** get the (imputed) missing values of a data set.
    *  @note In C++11, it will be possible to use a tuple rather that this
    *  pair of pair...
    *  @param data the array to return with the missing values
    **/
   template<typename Type_>
   void getMissingValues( std::vector< std::pair<std::pair<int,int>, Type_ > >& data) const;

   /** function to use in order to initialize super class and find
    *  missing values in the data set.
    **/
   void initialize()
   {
     // initialize dimensions
     this->setRows(dataij_.rows());
     this->setCols(dataij_.cols());
     nbVariable_ = this->sizeCols();
     // find coordinates of the missing values
     findMissing();
   }

   /** @return the element (i,j) of the 2D container.
    *  @param i, j indexes of the row and of the column
    **/
   inline Type& elt2Impl(int i, int j) { return dataij_.elt2Impl(i,j);}
   /** @return a constant reference of the element (i,j) of the 2D container.
    *  @param i, j indexes of the row and of the column
    **/
   inline ConstReturnType elt2Impl(int i, int j) const { return dataij_.elt2Impl(i,j);}
   /** @return a reference on the ith element
    *  @param i index of the ith element
    **/
   inline Type& elt1Impl(int i) { return dataij_.elt(i);}
   /** @return the constant ith element
    *  @param i index of the ith element
    **/
   inline ConstReturnType elt1Impl(int i) const { return dataij_.elt(i);}
   /** @return a reference on the number */
   inline Type& elt0Impl() { return dataij_.elt();}
   /** @return a constant reference on the number */
   inline ConstReturnType elt0Impl() const  { return dataij_.elt();}
   /** Access to the ith row of the Allocator.
    *  @param i index of the row
    *  @return a reference on the ith row
    **/
   inline Row row(int i) const { return dataij_.row(i);}
   /** Access to the row (i,J) of the Allocator.
    *  @param i,J index of the row and range of the columns
    *  @return a reference on the ith row
    **/
   inline SubRow row(int i, Range const& J) const { return dataij_.row(i, J);}
   /** Access to the jth column of the Allocator.
    *  @param j index of the column
    *  @return a reference on the jth column
    **/
   inline Col col(int j) const { return dataij_.col(j);}
   /** Access to the column (I,j) of the Allocator.
    *  @param I,j range of the rows and index of the column
    *  @return a reference on the jth column
    **/
   inline SubCol col(Range const& I, int j) const { return dataij_.col(I,j);}
   /** Access to the sub-part (I,J) of the Allocator.
    *  @param I,J range of the rows and columns
    *  @return a reference on a sub-part of the Allocator
    **/
   inline SubArray sub(Range const& I, Range const& J) const { return dataij_.sub(I,J);}
   /** Access to a sub-vector. For 1D allocators only.
    *  @param I range of the rows
    *  @return a reference on a sub-part of the Allocaor
    **/
   inline SubVector sub(Range const& I) const { return dataij_.sub(I);}

 protected:
   /** data set */
   Data dataij_;
   /** number of variables in the data set */
   int nbVariable_;
   /** utility function for lookup the data set and find missing values
    *  coordinates.
    *  @return the number of missing values
    **/
  virtual std::vector< std::pair<int,int> >::size_type findMissing();
};

template<class Data>
std::vector< std::pair<int,int> >::size_type DataBridge<Data>::findMissing()
{
#ifdef STK_DMANAGER_VERBOSE
  stk_cout << _T("IDataBridge::Entering findMissing()\n");
#endif
  for (int j=dataij_.beginCols(); j< dataij_.endCols(); ++j)
  {
    for (int i=dataij_.beginRows(); i< dataij_.endRows(); ++i)
    {
      if (Arithmetic<Type>::isNA(dataij_(i,j)))
      { v_missing_.push_back(std::pair<int,int>(i,j));}
    }
  }
 return v_missing_.size();
#ifdef STK_DMANAGER_VERBOSE
  stk_cout << _T("IDataBridge::findMissing() terminated, nbMiss= ") << v_missing_.size() << _T("\n");
#endif
}

template<class Data>
template<typename Type_>
void DataBridge<Data>::getMissingValues( std::vector< std::pair<std::pair<int,int>, Type_ > >& data) const
{
  data.resize(v_missing_.size());
  for(size_t i = 0; i< v_missing_.size(); ++i)
  {
    data[i].first  = v_missing_[i];
    data[i].second = dataij_(v_missing_[i].first, v_missing_[i].second);
  }
}
} // namespace STK

#endif /* STK_DATABRIDGE_H */
