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
 * Project:  rtkpp::
 * created on: 28 juil. 2014
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file RcppTraits.h
 *  @brief In this file .
 **/


#ifndef RCPPTRAITS_H
#define RCPPTRAITS_H

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 * Indicates the storage type associated with a SEXP type
 * for example: RcppTraits<int>::Rtype_ is a INTSXP
 * The default is VECSXP
 */
template<typename RTYPE> struct RcppTraits
{
  enum
  { Rtype_ = VECSXP};
};

/** Total specialization for integer vector (INTSXP)
 *  typedef from int (FIXME could be also logical)
 */
template<> struct RcppTraits<int>
{
  enum
  { Rtype_ = INTSXP};
};

/** Total specialization for numeric vectors (REALSXP)
 *  typedef from double
 */
template<> struct RcppTraits<double>
{
  enum
  { Rtype_ = REALSXP};
};

/** Total specialization for numeric vectors (CPLXSXP)
 *  typedef to Rcomplex
 */
template<> struct RcppTraits<Rcomplex>
{
  enum
  { Rtype_ = CPLXSXP};
};

/** Total specialization for raw vectors (RAWSXP)
 * typedef from Rbyte
 */
template<> struct RcppTraits<Rbyte>
{
  enum
  { Rtype_ = RAWSXP};
};

/*
 * Total specialization for logical vectors (LGLSXP)
 * typedef from int (!)
 */
//template<> struct RcppTraits<bool>
//{
//  enum
//  { Rtype_ = LGLSXP};
//};


} // namespace hidden

} // namespace STK

#endif /* RCPPTRAITS_H */
