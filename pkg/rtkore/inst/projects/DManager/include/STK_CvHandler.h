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
 * Project:  stkpp::DManager
 * created on: 15 nov. 2016
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_CvHandler.h
 *  @brief In this file we define the Interface Class for sampling data
 *  from data sets for CV and
 **/


#ifndef STK_CVHANDLER_H
#define STK_CVHANDLER_H

#include <Arrays/include/STK_CArrayVector.h>
#include <STatistiK/include/STK_Law_UniformDiscrete.h>

namespace STK
{

/** CvHanler is an utility function for building the submatrix/subvectors
 *  needed when using k-folds cross-validation.
 **/
class CvHandler
{
  public:
    /** Default constructor. nbFolds is set to the number of observation
     *  @param rangeData the range of the data to set
     *  @param nbFolds numbbe of Folds
     **/
    CvHandler( Range const& rangeData, int nbFolds)
             : rangeData_(rangeData), nbFolds_(nbFolds)
             , partitions_(rangeData), sizePartitions_(Range(0, nbFolds_), 0)
    { partition();}
    /** @return the number of folds */
    inline int nbFolds() const { return nbFolds_;}
    /** @return the range of the data */
    inline Range const& rangeData() const { return rangeData_;}
    /** @return the partitions */
    inline CVectorXi const& partitions() const { return partitions_;}
    /** @return the size of the partitions */
    inline CVectorXi const& sizePartitions() const { return sizePartitions_;}

    /** get the data set when setting out fold k and test data set  */
    template<class Data>
    void getKFold( int k, Data const& x,Data& xFold, Data& xTest) const
    {
      // check dimensions
      if (x.rows() != rangeData_)
      { STKRUNTIME_ERROR_1ARG(CvHandler::getKFold,k,x.rows()!=rangeData_);}
      if (sizePartitions_.begin() > k)
      { STKOUT_OF_RANGE_1ARG(CvHandler::getKFold,k,k<sizePartitions_.begin());}
      if (sizePartitions_.end() <= k)
      { STKOUT_OF_RANGE_1ARG(CvHandler::getKFold,k,k<sizePartitions_.end()<=k);}
      // prepare containers
      Range xFoldRows = x.rows();
      xFoldRows.decLast(sizePartitions_[k]);
      xFold.resize(xFoldRows, x.cols());
      xTest.resize(sizePartitions_[k], x.cols());
      // copy data
      int iFoldRow = xFold.beginRows(), iTestRow = xTest.beginRows();
      for (int i = partitions_.begin(); i < partitions_.end(); ++i)
      {
        if (partitions_[i] == k)
        {
          xTest.row(iTestRow) = x.row(i);
          ++iTestRow;
        }
        else
        {
          xFold.row(iFoldRow) = x.row(i);
          ++iFoldRow;
        }
      }
    }
    /** get the data set when setting out fold k and test data set  */
    template<class xData, class yData>
    void getKFold( int k, xData const& x, xData& xFold, xData& xTest
                        , yData const& y, yData& yFold, yData& yTest) const
    {
      // check dimensions
      if (x.rows() != rangeData_)
      { STKRUNTIME_ERROR_1ARG(CvHandler::getKFold,k,x.rows()!=rangeData_);}
      // check dimensions
      if (y.rows() != rangeData_)
      { STKRUNTIME_ERROR_1ARG(CvHandler::getKFold,k,y.rows()!=rangeData_);}
      if (sizePartitions_.begin() > k)
      { STKOUT_OF_RANGE_1ARG(CvHandler::getKFold,k,k<sizePartitions_.begin());}
      if (sizePartitions_.end() <= k)
      { STKOUT_OF_RANGE_1ARG(CvHandler::getKFold,k,k<sizePartitions_.end()<=k);}
      // prepare constainers
      Range xFoldRows = x.rows();
      xFoldRows.decLast(sizePartitions_[k]);
      xFold.resize(xFoldRows, x.cols());
      xTest.resize(sizePartitions_[k], x.cols());
      yFold.resize(xFoldRows, y.cols());
      yTest.resize(sizePartitions_[k], y.cols());
      // copy data
      int iFoldRow = xFold.beginRows(), iTestRow = xTest.beginRows();
      for (int i = partitions_.begin(); i < partitions_.end(); ++i)
      {
        if (partitions_[i] == k)
        {
          xTest.row(iTestRow) = x.row(i);
          yTest.row(iTestRow) = y.row(i);
          ++iTestRow;
        }
        else
        {
          xFold.row(iFoldRow) = x.row(i);
          yFold.row(iFoldRow) = y.row(i);
          ++iFoldRow;
        }
      }
    }
  protected:
    /** create a random partition in k folds*/
    inline void partition()
    {
      //fill the container with the index of folds
      for(int i = partitions_.begin() ; i< partitions_.end() ;i++)
      {
        partitions_[i] = i%nbFolds_;
        sizePartitions_[i%nbFolds_]++;
      }
      //make a random rearrangement
      int begin = partitions_.begin();
      for (int i=partitions_.end()-2; i>begin; --i)
      { std::swap(partitions_[i], partitions_[Law::UniformDiscrete::rand(begin, i+1)]);}
    }

  private:
    /** Range of the data set (number of rows) */
    Range rangeData_;
    /** Number of folds */
    int nbFolds_;
    /** repartition of the sample into k-folds */
    CVectorXi partitions_;
    /** size of each fold */
    CVectorXi sizePartitions_;
};

} // namespace STK

#endif /* STK_CVHANDLER_H */
