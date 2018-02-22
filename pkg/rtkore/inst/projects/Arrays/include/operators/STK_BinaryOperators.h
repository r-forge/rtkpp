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
 * Project:  stkpp::Arrays
 * created on: 17 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_BinaryOperators.h
 *  @brief In this file we implement the BinaryOperator class.
 **/


#ifndef STK_BINARYOPERATORS_H
#define STK_BINARYOPERATORS_H

#include "STK_SlicingOperators.h"

#define EGAL(arg1, arg2) ((arg1::structure_ == int(Arrays::arg2)))

namespace STK
{

namespace hidden
{
/** @ingroup hidden
 *  @brief Helper class giving the Structure of a binary operator
 **/
template<typename FunctorOp, typename Lhs, typename Rhs, int LStructure_, int RStructure_>
struct BinaryEltImpl;

//----------------
// Lhs array2D_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_      = Arrays::array2D_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::array2D_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_2D_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};



//----------------
// Lhs square_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_ = Arrays::binary_op_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::square_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_2D_LowSym_
       , sizeRows_  = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_  = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};


//----------------
// lhs diagonal_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_Diag_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Diag_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_Diag_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f( LType(), r.elt(i,j)) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_Diag_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f( LType(), r.elt(i,j)) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Diag_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::upper_symmetric_
       , binary_op_Kind_= Arrays::binary_op_Diag_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f( LType(), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::lower_symmetric_
       , binary_op_Kind_= Arrays::binary_op_Diag_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f( LType(), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};



//-------------------------
// lhs upper_triangular_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_UpTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_UpTri_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_UpTri_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_UpTri_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_triangular_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpTri_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(j,i)) : f( LType(), r.elt(i,j));}
};



// lhs is lower_triangular_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_LowTri_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::array2D_
       , binary_op_Kind_= Arrays::binary_op_LowTri_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_LowTri_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(j,i)) : f( LType(), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_triangular_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( LType(), r.elt(j,i));}
};



// lhs is symmetric_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_Sym_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::symmetric_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_Sym_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(i,j), r.elt(j,i));}
};


//---------------------
// lhs upper_symmetric_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpSym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpSym_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
      };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::upper_triangular_
       , binary_op_Kind_= Arrays::binary_op_UpSym_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_UpSym_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i<j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::upper_symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j,i), r.elt(j,i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::upper_symmetric_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_UpSym_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i<=j) ? f(l.elt(i,j), r.elt(j,i)) : f( l.elt(j,i), r.elt(i,j));}
};


//-----------------------
// lhs lower_symmetric_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::array2D_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j, i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::square_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowTri_2D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j, i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowTri_Diag_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::upper_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::square_
       , binary_op_Kind_= Arrays::binary_op_LowSym_UpTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : (i>j) ? f(l.elt(i,j), RType()) : f(l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::lower_triangular_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::lower_triangular_
       , binary_op_Kind_= Arrays::binary_op_LowSym_LowTri_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(j,i), RType());}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowSym_Sym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f( l.elt(j, i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::upper_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowSym_UpSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(j,i)) : f(l.elt(j,i), r.elt(i,j));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::lower_symmetric_, Arrays::lower_symmetric_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::lower_symmetric_
       , binary_op_Kind_= Arrays::binary_op_LowSym_LowSym_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i>=j) ? f(l.elt(i,j), r.elt(i,j)) : f(l.elt(j,i), r.elt(j,i));}
};


//--------------------------------------------------------------------------------------------
// 1D case

//-----------------------
// Lhs diagonal_, Rhs 1D
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  typedef typename Lhs::Type LType;
  typedef typename Rhs::Type RType;
  enum { structure_ = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return (i==j) ? f(l.elt(i,j), r.elt(i,j)) : f(LType(), RType());}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::point_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeCols_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::diagonal_, Arrays::vector_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_   = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_   = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeRows_)
       , useLhsRows_ = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_ = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};

//-----------------------
// Lhs vector_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::vector_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Rhs::sizeRows_ != UnknownSize ? int(Rhs::sizeRows_) : int(Lhs::sizeCols_)
       , sizeCols_      = Rhs::sizeCols_ != UnknownSize ? int(Rhs::sizeCols_) : int(Lhs::sizeCols_)
       , useLhsRows_    = Rhs::sizeRows_ != UnknownSize ? Arrays::useRhsSize_ : Arrays::useLhsSize_
       , useLhsCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::vector_, Arrays::vector_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::vector_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeRows_)
       , sizeCols_      = 1
       , useLhsRows_    = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       , useLhsCols_    = Arrays::useLhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::vector_, Arrays::point_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::vector_ // vector + point is vector (somehow an arbitrary choice)
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Lhs::sizeRows_ != UnknownSize ? int(Lhs::sizeRows_) : int(Rhs::sizeCols_)
       , sizeCols_      = 1
       , useLhsRows_    = Lhs::sizeRows_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       , useLhsCols_    = Arrays::useLhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};

//-----------------------
// Lhs point_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::point_, Arrays::diagonal_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::diagonal_ // point_ + diagonal_ is diagonal_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = Rhs::sizeRows_ != UnknownSize ? int(Rhs::sizeRows_) : int(Lhs::sizeCols_)
       , sizeCols_      = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_    = Rhs::sizeRows_ != UnknownSize ? Arrays::useRhsSize_ : Arrays::useLhsOtherSize_
       , useLhsCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::point_, Arrays::vector_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::point_ // point + vector_ is point_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = 1
       , sizeCols_      = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeRows_)
       , useLhsRows_    = Arrays::useLhsSize_
       , useLhsCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsOtherSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::point_, Arrays::point_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_ = Arrays::point_
       , binary_op_Kind_= Arrays::binary_op_1D_
       , sizeRows_      = 1
       , sizeCols_      = Lhs::sizeCols_ != UnknownSize ? int(Lhs::sizeCols_) : int(Rhs::sizeCols_)
       , useLhsRows_    = Arrays::useLhsSize_
       , useLhsCols_    = Lhs::sizeCols_ != UnknownSize ? Arrays::useLhsSize_ : Arrays::useRhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
};

//-----------------------
// Lhs number_
template< typename FunctorOp, typename Lhs, typename Rhs>
struct BinaryEltImpl< FunctorOp, Lhs, Rhs, Arrays::number_, Arrays::number_>
{
  typedef typename FunctorOp::result_type result_type;
  enum { structure_     = Arrays::number_
       , binary_op_Kind_= Arrays::binary_op_0D_
       , sizeRows_      = 1
       , sizeCols_      = 1
       , useLhsRows_    = Arrays::useLhsSize_
       , useLhsCols_    = Arrays::useLhsSize_
       };
  inline static result_type elt2Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i, int j)
  { return f(l.elt(i,j), r.elt(i,j));}
  inline static result_type elt1Impl(FunctorOp const& f, Lhs const& l, Rhs const& r, int i)
  { return f(l.elt(i), r.elt(i));}
  inline static result_type elt0Impl(FunctorOp const& f, Lhs const& l, Rhs const& r)
  { return f(l.elt(), r.elt());}
};

} //namespace hidden


// forward declaration
template<typename FunctorOp, typename Lhs, typename Rhs>
class BinaryOperator;

namespace hidden {

/** @ingroup hidden
 *  @brief implement the access to the rows of the BinaryOperator
 *  Possible cases are:
 *  - use lhs.rows()
 *  - use rhs.rows()
 *  - use lhs.cols()
 *  - use rhs.cols()
 **/
template< typename Lhs, typename Rhs, int Size_, int useLhsRows_>
struct BinaryRowsImpl;
/** @ingroup hidden
 *  @brief implement the access to the columns of the BinaryOperator
 *  Possible cases are:
 *  - use lhs.cols()
 *  - use rhs.cols()
 *  - use lhs.rows()
 **/
template< typename Lhs, typename Rhs, int Size_, int useLhsCols_>
struct BinaryColsImpl;
/** @ingroup hidden
  * @brief specialization for the case useLhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl< Lhs, Rhs,  Size_, Arrays::useLhsSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.rows();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl<Lhs, Rhs, Size_, Arrays::useRhsSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.rows();}
};
/** @ingroup hidden
  * @brief specialization for the case useLhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl<Lhs, Rhs, Size_, Arrays::useLhsOtherSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryRowsImpl<Lhs, Rhs, Size_, Arrays::useRhsOtherSize_>
{
  /** Type of the Range for the rows */
  typedef TRange<Size_> RowRange;
  /**  @return the range of the rows */
  inline static RowRange const& rowsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useLhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl<Lhs, Rhs, Size_, Arrays::useLhsSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl< Lhs, Rhs, Size_, Arrays::useRhsSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.cols();}
};
/** @ingroup hidden
  * @brief specialization for the case useLhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl< Lhs, Rhs, Size_, Arrays::useLhsOtherSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return lhs.rows();}
};
/** @ingroup hidden
  * @brief specialization for the case useRhsOtherSize_
  **/
template<typename Lhs, typename Rhs, int Size_>
struct BinaryColsImpl< Lhs, Rhs,Size_, Arrays::useRhsOtherSize_>
{
  /** Type of the Range for the columns */
  typedef TRange<Size_> ColRange;
  /**  @return the range of the columns */
  inline static ColRange const& colsImpl(Lhs const& lhs, Rhs const& rhs) { return rhs.rows();}
};


/** @ingroup hidden
 *  @brief Traits class for the BinaryOperator
 */
template<typename FunctorOp, typename Lhs, typename Rhs>
struct Traits< BinaryOperator<FunctorOp, Lhs, Rhs> >
{
  enum
  {
    // find the kind of binary operator and the Structure using helper class BinaryEltImpl
    binary_op_Kind_ = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::binary_op_Kind_,

    isLhs1D_ = EGAL(Lhs,vector_)||EGAL(Lhs,point_)||EGAL(Lhs,diagonal_),
    isRhs1D_ = EGAL(Rhs,vector_)||EGAL(Rhs,point_)||EGAL(Rhs,diagonal_),

    isRhs2D_ = EGAL(Rhs,array2D_)||EGAL(Rhs,square_)||EGAL(Rhs,diagonal_)
             ||EGAL(Rhs,lower_triangular_)||EGAL(Rhs,upper_triangular_)
             ||EGAL(Rhs,symmetric_)||EGAL(Rhs,lower_symmetric_)||EGAL(Rhs,upper_symmetric_),
    isLhs2D_ = EGAL(Lhs,array2D_)||EGAL(Lhs,square_)||EGAL(Lhs,diagonal_)
             ||EGAL(Lhs,lower_triangular_)||EGAL(Lhs,upper_triangular_)
             ||EGAL(Lhs,symmetric_)||EGAL(Lhs,lower_symmetric_)||EGAL(Lhs,upper_symmetric_),


    isRes0D_ = EGAL(Lhs,number_) && EGAL(Rhs,number_),
    isRes1D_ = (EGAL(Lhs,vector_)||EGAL(Lhs,point_)) && (EGAL(Rhs,vector_)||EGAL(Rhs,point_)),
    isRes2D_ = isLhs2D_ && isRhs2D_,

    is1D1D_  = isLhs1D_ && isRhs1D_,

    // get the structure from the helper class BinaryEltImpl
    structure_ = hidden::BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::structure_,
    orient_    = Lhs::orient_,    // preserve the Lhs storage orientation. Could be optimized ?
    sizeRows_  = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::sizeRows_,
    sizeCols_  = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::sizeCols_,
    storage_   = (Lhs::storage_ == int(Arrays::dense_)) || (Rhs::storage_ == int(Arrays::dense_))
               ?  int(Arrays::dense_) : int(Arrays::sparse_),

    useLhsRows_    = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::useLhsRows_,
    useLhsCols_    = BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_>::useLhsCols_
  };
  typedef RowOperator<BinaryOperator<FunctorOp, Lhs, Rhs> > Row;
  typedef ColOperator<BinaryOperator<FunctorOp, Lhs, Rhs> > Col;
  typedef typename FunctorOp::result_type Type;
  typedef typename FunctorOp::result_type ConstReturnType;
};

} // end namespace hidden


/** @ingroup Arrays
  * @brief Generic expression where a binary operator is applied to two expressions
  *
  * @tparam FunctorOp template functor implementing the binary operator
  * @tparam Lhs left-hand side type
  * @tparam Rhs right-hand side type
  *
  * This class represents an expression  where a binary operator is applied to
  * two expressions.
  * It is the return type of binary operators, by which we mean only those
  * binary operators where both the left-hand side and the right-hand side
  * are expressions. For example, the return type of matrix1+matrix2 is a
  * BinaryOperator. The return type of number+matrix is a unary operator.
  *
  * Most of the time, this is the only way that it is used, so you typically
  * don't have to name BinaryOperator types explicitly.
  **/

template<typename FunctorOp, typename Lhs, typename Rhs>
class BinaryOperator: public ExprBase< BinaryOperator<FunctorOp, Lhs, Rhs> >
                    , public TRef<1>
{
  public:
    typedef hidden::BinaryEltImpl<FunctorOp, Lhs, Rhs, Lhs::structure_, Rhs::structure_> EltImpl;
    typedef hidden::BinaryRowsImpl< Lhs, Rhs, EltImpl::sizeRows_, EltImpl::useLhsRows_ > RowsImpl;
    typedef hidden::BinaryColsImpl< Lhs, Rhs, EltImpl::sizeCols_, EltImpl::useLhsCols_ > ColsImpl;

    typedef typename hidden::Traits< BinaryOperator >::ConstReturnType ConstReturnType;
    typedef typename hidden::Traits<BinaryOperator >::Type Type;
//    typedef typename hidden::Traits<BinaryOperator >::Row Row;
//    typedef typename hidden::Traits<BinaryOperator >::Col Col;

    enum
    {
      isRes0D_   = hidden::Traits<BinaryOperator>::isRes0D_,
      isRes1D_   = hidden::Traits<BinaryOperator>::isRes1D_,
      isRes2D_   = hidden::Traits<BinaryOperator>::isRes2D_,
      is1D1D_    = hidden::Traits<BinaryOperator>::is1D1D_,

      structure_ = hidden::Traits<BinaryOperator>::structure_,
      orient_    = hidden::Traits<BinaryOperator>::orient_,
      sizeRows_  = hidden::Traits<BinaryOperator>::sizeRows_,
      sizeCols_  = hidden::Traits<BinaryOperator>::sizeCols_,
      storage_   = hidden::Traits<BinaryOperator>::storage_,

      // All the valid cases for binary operators
      isValid_ =( isRes0D_ || isRes1D_  || isRes2D_ || is1D1D_)
    };
    /** Type of the Range for the rows */
    typedef TRange<sizeRows_> RowRange;
    /** Type of the Range for the columns */
    typedef TRange<sizeCols_> ColRange;

    /** default constructor */

    inline BinaryOperator( Lhs const& lhs, Rhs const& rhs, FunctorOp const& func = FunctorOp())
                         : lhs_(lhs), rhs_(rhs), functor_(func)
    {
      // FIXME : not safe. Add more test in the 1D case at compile time (and runtime ?)
      STK_STATIC_ASSERT_BINARY_OPERATOR_MISMATCH( isValid_ );
      STK_STATIC_ASSERT_COLS_DIMENSIONS_MISMATCH(!( (int(Lhs::sizeCols_) != UnknownSize)
                                                &&  (int(Rhs::sizeCols_) != UnknownSize)
                                                &&  (int(Lhs::sizeCols_) != int(Rhs::sizeCols_))
                                                &&  (isRes2D_)
                                                 ));
      STK_STATIC_ASSERT_ROWS_DIMENSIONS_MISMATCH(!( (int(Lhs::sizeRows_) != UnknownSize)
                                                &&  (int(Rhs::sizeRows_) != UnknownSize)
                                                &&  (int(Lhs::sizeRows_) != int(Rhs::sizeRows_))
                                                &&  (isRes2D_)
                                                 ));
#ifdef STK_BOUNDS_CHECK
      if ((lhs.rows() != rhs.rows()) && (isRes2D_))
      { STKRUNTIME_ERROR_2ARG(BinaryOperator, lhs.rows(), rhs.rows(), Rows range mismatch in BinaryOperator);}
      if (( lhs.cols() != rhs.cols()) && (isRes2D_))
      { STKRUNTIME_ERROR_2ARG(BinaryOperator, lhs.cols(), rhs.cols(), Columns range mismatch in BinaryOperator);}
#endif
    }

    /** @return the left hand side expression */
    inline Lhs const& lhs() const { return lhs_; }
    /** @return the right hand side nested expression */
    inline Rhs const& rhs() const { return rhs_; }
    /** @return the functor representing the binary operation */
    inline FunctorOp const& functor() const { return functor_; }

    /** @return element (i,j) */
    inline ConstReturnType elt2Impl(int i, int j) const { return EltImpl::elt2Impl(functor_, lhs_, rhs_, i, j);}
    /** @return element i */
    inline ConstReturnType elt1Impl(int i) const { return EltImpl::elt1Impl(functor_, lhs_, rhs_, i);}
    /** @return element */
    inline ConstReturnType elt0Impl() const { return EltImpl::elt0Impl(functor_, lhs_, rhs_);}
    /** @return range of the rows */
    inline RowRange const& rowsImpl() const { return RowsImpl::rowsImpl(lhs_, rhs_);}
    /** @return range of the columns */
    inline ColRange const& colsImpl() const { return ColsImpl::colsImpl(lhs_, rhs_);}

  protected:
    Lhs const& lhs_;
    Rhs const& rhs_;
    FunctorOp const functor_;
};


} // namespace STK

#endif /* STK_BINARYOPERATORS_H */
