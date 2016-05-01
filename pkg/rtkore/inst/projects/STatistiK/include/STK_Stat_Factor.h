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
 * Project:  stkpp::STatistiK
 * Purpose:  Compute elementary 1D statistics for all variables.
 * Author:   Serge Iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_Stat_Factor.h
 *  @brief In this file we define and implement the Factor class.
 **/

#ifndef STK_STAT_FACTOR_H
#define STK_STAT_FACTOR_H

#include <Sdk/include/STK_IRunner.h>

namespace STK
{
namespace Stat
{
/** @ingroup StatDesc
 *  @brief Computation of the Factors of a 2D Container.
 *
 *  The class @c Factor is a factory class for computing the factors of an array.
 *  The values can be of any type. the Coding is performed from the previous type
 *  in integer. The mapping is stored and preserved in an array of map array.
 **/
template <class Array>
class Factor: public IRunnerWithData<Array>
{
  public:
    typedef IRunnerWithData<Array> Base;
    typedef typename Array::Row RowVector;
    typedef typename Array::Col ColVector;
    typedef typename Array::Type Type;

    typedef std::map<Type, int> EncodeMap;
    typedef CArrayPoint<EncodeMap> Encoder;
    using Base::p_data_;

    /** Default Constructor. */
    Factor(): Base(), levels_(), nbLevels_(), coding_()
    {}
    /** Constructor.
     *  @param data a reference on the data set
     **/
    Factor( Array const& data): Base(data)
                              , levels_(p_data_->rows(),p_data_->cols())
                              , nbLevels_(p_data_->cols(), baseIdx - 1)
                              , coding_(p_data_->cols())

    {}
    /** Constructor.
     *  @param p_data a pointer on the data set
     **/
    Factor( Array const* p_data): Base(p_data), levels_(), nbLevels_(), coding_()
    {
      if (p_data_)
      {
        levels_.resize(p_data_->rows(),p_data_->cols());
        nbLevels_.resize(p_data_->cols()).setValue(baseIdx - 1);
        coding_.resize(p_data_->cols());
      }
    }
    /** copy constructor.
     *  @param f the Factor to copy
     **/
    Factor( Factor const& f): Base(f), levels_(f.levels_), nbLevels_(f.nbLevels_)
                            , coding_(f.coding_)
    {}
    /** virtual destructor.*/
    virtual ~Factor() {}

    /** @return the array of factor levels */
    inline CArrayXXi const& levels() const { return levels_;}
    /** @return the array of number of levels */
    inline CPointXi const& nbLevels() const { return nbLevels_;}
    /** @return the array of encoding map */
    inline Encoder const& coding() const { return coding_;}

    /** clone pattern */
    inline virtual Factor* clone() const { return new Factor(*this);}
    /** run the estimation of the Factor statistics. **/
    virtual bool run();

  protected:
    /** array of levels */
    CArrayXXi levels_;
    /** Number of levels of each variables */
    CPointXi nbLevels_;
    /** coding of the original variables */
    Encoder coding_;

    /** udpating method in case we set a new data set */
    virtual void update();
};

template <class Array>
void Factor<Array>::update()
{
   // if there is no data there is nothing to update
   if (p_data_)
   {
     levels_.resize(p_data_->rows(),p_data_->cols());
     nbLevels_.resize(p_data_->cols()).setValue(baseIdx - 1);
     for(int j=coding_.begin(); j<std::min(coding_.end(), p_data_->cols().end()); ++j)
     { coding_[j].clear();}
     coding_.resize(p_data_->cols());
   }
}

template <class Array>
bool Factor<Array>::run()
{
  if (!p_data_)
  { this->msg_error_ = STKERROR_NO_ARG(FactorArray::run,data is not set);
    return false;
  }
  try
  {
    for (int j=p_data_->beginCols(); j< p_data_->endCols(); ++j)
    {
      for (int i=p_data_->beginRows(); i< p_data_->endRows(); ++i)
      {
        // find coding
        Type idData = p_data_->elt(i,j);
        typename EncodeMap::const_iterator it = coding_[j].find(idData);
        if (it != coding_[j].end()) { levels_(i,j) = it->second;}
        else // find a new level to add
        {
          levels_(i,j) = (++nbLevels_[j]); // create a new level and set it
          coding_[j].insert(std::pair<Type, int>(idData, nbLevels_[j]));
        }
      }
    }
  }
  catch (Exception const& error)
  {
    this->msg_error_ += _T("Error in Factor::run():\nWhat: ");
    this->msg_error_ += error.error();
    return false;
  }
  // no error
  return true;
}

}  // namespace Stat

}  // namespace STK

#endif /*STK_STAT_FACTOR_H */
