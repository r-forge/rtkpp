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
 * created on: 16 oct. 2012
 * Author:   iovleff, S..._Dot_I..._At_stkpp_Dot_org (see copyright for ...)
 **/

/** @file STK_AssignSelector.h
 *  @brief In this file we implement the copy and assign methods used when
 *  copying an array or an expression (rhs) in an other array (lhs).
 **/


#ifndef STK_ASSIGNSELECTOR_H
#define STK_ASSIGNSELECTOR_H


namespace STK
{


namespace hidden
{
/** @ingroup hidden
 *  @brief Copycat to use at compile time.
 *  If n is the number of structures, there is potentially n^2 ways to
 *  copy rhs inside lhs.
 *  */
template < typename Derived, typename Rhs, int TStructure_, int RhsStructure_>
struct Copycat;

//---------------------GENERAL----------------------------------
// general <- general
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::array2D_, Arrays::array2D_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    for (int j = rhs.beginCols(); j < rhs.endCols(); ++j)
      for (int i = rhs.beginRows(); i < rhs.endRows(); ++i)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    for (int i = rhs.beginRows(); i < rhs.endRows(); ++i)
      for (int j = rhs.beginCols(); j < rhs.endCols(); ++j)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
  static void copyByCol(Lhs& lhs, Rhs const& rhs )
  {
    for ( int jRhs=rhs.beginCols(), jLhs=lhs.beginCols(); jRhs<rhs.endCols(); jLhs++, jRhs++)
      for ( int iRhs=rhs.beginRows(), iLhs=lhs.beginRows(); iRhs<rhs.endRows(); iLhs++, iRhs++)
    { lhs.elt(iLhs, jLhs) = rhs.elt(iRhs, jRhs);}
  }
  static void copyByRow(Lhs& lhs, Rhs const& rhs )
  {
    for ( int iRhs=rhs.beginRows(), iLhs=lhs.beginRows(); iRhs<rhs.endRows(); iLhs++, iRhs++)
      for ( int jRhs=rhs.beginCols(), jLhs=lhs.beginCols(); jRhs<rhs.endCols(); jLhs++, jRhs++)
      { lhs.elt(iLhs, jLhs) = rhs.elt(iRhs, jRhs);}
  }
};

// general <- square
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::array2D_, Arrays::square_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    for (int i = rhs.begin(); i< rhs.end(); ++i)
      for (int j = rhs.begin(); j < rhs.end(); ++j)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    for (int j = rhs.begin(); j < rhs.lend(); ++j)
      for (int i = rhs.begin(); i< rhs.end(); ++i)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
  static void copyByCol(Lhs& lhs, Rhs const& rhs )
  {
    for ( int jRhs=rhs.beginCols(), jLhs=lhs.beginCols(); jRhs<rhs.endCols(); jLhs++, jRhs++)
      for ( int iRhs=rhs.beginRows(), iLhs=lhs.beginRows(); iRhs<rhs.endRows(); iLhs++, iRhs++)
    { lhs.elt(iLhs, jLhs) = rhs.elt(iRhs, jRhs);}
  }
  static void copyByRow(Lhs& lhs, Rhs const& rhs )
  {
    for ( int iRhs=rhs.beginRows(), iLhs=lhs.beginRows(); iRhs<rhs.endRows(); iLhs++, iRhs++)
      for ( int jRhs=rhs.beginCols(), jLhs=lhs.beginCols(); jRhs<rhs.endCols(); jLhs++, jRhs++)
      { lhs.elt(iLhs, jLhs) = rhs.elt(iRhs, jRhs);}
  }
};

// general <- diagonal
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::array2D_, Arrays::diagonal_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    lhs.setValue(Type(0));
    for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i,i) = rhs.elt(i);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    lhs.setValue(Type(0));
    for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i,i) = rhs.elt(i);}
  }
};

// general <- lower_triangular
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::array2D_,   Arrays::lower_triangular_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int j = rhs.beginCols(); j < end; ++j)
    { int i = rhs.beginRows();
      for (; i < j; ++i)             { lhs.elt(i,j) = Type(0);}
      for (; i < rhs.endRows(); ++i) { lhs.elt(i,j) = rhs.elt(i, j);}
    }
    for (int j= end; j < rhs.endCols(); ++j)
      for (int i = rhs.beginRows(); i < rhs.endRows(); ++i) { lhs.elt(i,j) = Type(0);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int i = rhs.beginRows(); i < end; ++i)
    { int j = rhs.beginCols();
      for (; j <=i; ++j)             { lhs.elt(i,j) = rhs.elt(i, j);}
      for (; j < rhs.endCols(); ++j) { lhs.elt(i,j) = Type(0);}
    }
    for (int i= end; i < rhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j < rhs.endCols(); ++j) { lhs.elt(i,j) = rhs.elt(i, j);}
  }
};

// general <- upper_triangular
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::array2D_, Arrays::upper_triangular_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int j = rhs.beginCols(); j < end; ++j)
    { int i = rhs.beginRows();
      for (; i <= j; ++i)            { lhs.elt(i,j) = rhs.elt(i, j);}
      for (; i < rhs.endRows(); ++i) { lhs.elt(i,j) =  Type(0);}
    }
    for (int j= end; j < rhs.endCols(); ++j)
      for (int i = rhs.beginRows(); i < rhs.endRows(); ++i) { lhs.elt(i,j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(lhs.endRows(), lhs.endCols());
    for (int i = rhs.beginRows(); i < end; ++i)
    { int j = rhs.beginCols();
      for (; j <i; ++j)               { lhs.elt(i,j) = Type(0);}
      for (; j < rhs.endCols(); ++j) { lhs.elt(i,j) = rhs.elt(i, j);}
    }
    for (int i= end; i < rhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j < rhs.endCols(); ++j) { lhs.elt(i,j) = Type(0);}
  }
};

// general <- vector
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::array2D_, Arrays::vector_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    int j = rhs.beginCols();
    for (int i = rhs.beginRows(); i < rhs.endRows(); ++i)
    { lhs.elt(i, j) = rhs.elt(i);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    int j = rhs.beginCols();
    for (int i = rhs.beginRows(); i < rhs.endRows(); ++i)
    { lhs.elt(i, j) = rhs.elt(i);}
  }
};

// general <- point
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::array2D_, Arrays::point_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    int i = rhs.beginRows();
    for (int j = rhs.beginCols(); j < rhs.endCols(); ++j)
    { lhs.elt(i, j) = rhs.elt(j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    int i = rhs.beginRows();
    for (int j = rhs.beginCols(); j < rhs.endCols(); ++j)
    { lhs.elt(i, j) = rhs.elt(j);}
  }
};

//---------------------SQUARE----------------------------------
// square <- general
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::square_, Arrays::array2D_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    for (int i = rhs.beginRows(); i < rhs.endRows(); ++i)
      for (int j = rhs.beginCols(); j < rhs.endCols(); ++j)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    for (int j = rhs.beginCols(); j < rhs.endCols(); ++j)
      for (int i = rhs.beginRows(); i < rhs.endRows(); ++i)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
};

// square <- square
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::square_, Arrays::square_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    for (int j = rhs.begin(); j < rhs.end(); ++j)
      for (int i = rhs.begin(); i< rhs.end(); ++i)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    for (int i = rhs.begin(); i< rhs.end(); ++i)
      for (int j = rhs.begin(); j < rhs.end(); ++j)
      { lhs.elt(i, j) = rhs.elt(i, j);}
  }
};

// square <- diagonal
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::square_, Arrays::diagonal_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    lhs.setValue(Type(0));
    const int end = std::min(lhs.endRows(), lhs.endCols());
    for (int i = rhs.beginRows(); i < end; ++i) { lhs.elt(i,i) = rhs.elt(i);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    lhs.setValue(Type(0));
    const int end = std::min(lhs.endRows(), lhs.endCols());
    for (int i = rhs.beginRows(); i < end; ++i) { lhs.elt(i,i) = rhs.elt(i);}
  }
};

// square_ <- lower_triangular
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::square_,   Arrays::lower_triangular_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int j = rhs.beginCols(); j < end; ++j)
    { int i = rhs.beginRows();
      for (; i < j; ++i)                 { lhs.elt(i,j) = Type(0);}
      for (; i < rhs.endRows(); ++i) { lhs.elt(i,j) = rhs.elt(i, j);}
    }
    for (int j= end; j < rhs.endCols(); ++j)
      for (int i = rhs.beginRows(); i < rhs.endRows(); ++i) { lhs.elt(i,j) = Type(0);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int i = rhs.beginRows(); i < end; ++i)
    { int j = rhs.beginCols();
      for (; j <=i; ++j)                 { lhs.elt(i,j) = rhs.elt(i, j);}
      for (; j < rhs.endCols(); ++j) { lhs.elt(i,j) = Type(0);}
    }
    for (int i= end; i < rhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j < rhs.endCols(); ++j) { lhs.elt(i,j) = rhs.elt(i, j);}
  }
};

// square_ <- upper triangular
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::square_, Arrays::upper_triangular_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int j = rhs.beginCols(); j < end; ++j)
    { int i = rhs.beginRows();
      for (; i <= j; ++i)            { lhs.elt(i,j) = rhs.elt(i, j);}
      for (; i < rhs.endRows(); ++i) { lhs.elt(i,j) =  Type(0);}
    }
    for (int j= end; j < rhs.endCols(); ++j)
      for (int i = rhs.beginRows(); i < rhs.endRows(); ++i) { lhs.elt(i,j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(lhs.endRows(), lhs.endCols());
    for (int i = rhs.beginRows(); i < end; ++i)
    { int j = rhs.beginCols();
      for (; j <i; ++j)                  { lhs.elt(i,j) = Type(0);}
      for (; j < rhs.endCols(); ++j) { lhs.elt(i,j) = rhs.elt(i, j);}
    }
    for (int i= end; i < rhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j < rhs.endCols(); ++j) { lhs.elt(i,j) = Type(0);}
  }
};

//---------------------LDO----------------------------------
// lower_triangular <- lower_triangular
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::lower_triangular_, Arrays::lower_triangular_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int j = rhs.beginCols(); j < end; ++j)
      for (int i=j; i < rhs.endRows(); ++i)
      { lhs.elt(i,j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int i = rhs.beginRows(); i < end; ++i)
    {
      for (int j = rhs.beginCols(); j <=i; ++j) { lhs.elt(i,j) = rhs.elt(i, j);}
    }
    for (int i= end; i < rhs.endRows(); ++i)
      for (int j=rhs.beginCols(); j < rhs.endCols(); ++j) { lhs.elt(i,j) = rhs.elt(i, j);}
  }
};

//---------------------LUP----------------------------------
// upper_triangular <- upper_triangular
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::upper_triangular_, Arrays::upper_triangular_>
{
  typedef typename hidden::Traits<Lhs>::Type Type;
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  {
    const int end = std::min(rhs.endRows(), rhs.endCols());
    for (int j = rhs.beginCols(); j < end; ++j)
    {
      for (int i = rhs.beginRows(); i <= j; ++i) { lhs.elt(i,j) = rhs.elt(i, j);}
    }
    for (int j= end; j < rhs.endCols(); ++j)
      for (int i = rhs.beginRows(); i < rhs.endRows(); ++i) { lhs.elt(i,j) = rhs.elt(i, j);}
  }
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  {
    const int last = std::min(lhs.lastIdxRows(), lhs.lastIdxCols());
    for (int i = rhs.beginRows(); i <= last; ++i)
    { for (int j=i; j < rhs.endCols(); ++j) { lhs.elt(i,j) = rhs.elt(i, j);}}
  }
};


//---------------------DIAGONAL----------------------------------
// diagonal <- diagonal
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::diagonal_, Arrays::diagonal_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};

// diagonal <- vector
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::diagonal_, Arrays::vector_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};

// diagonal <- point
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::diagonal_, Arrays::point_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};

//---------------------VECTOR----------------------------------
//  vector <- diagonal
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::vector_, Arrays::diagonal_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};

//  vector <- vector
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::vector_, Arrays::vector_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};
// vector <- point
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::vector_, Arrays::point_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};

//---------------------POINT----------------------------------
//  point_ <- diagonal
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::point_, Arrays::diagonal_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};

//  vector <- vector
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::point_, Arrays::vector_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};
// vector <- point
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::point_, Arrays::point_>
{
  static void runByCol(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
  static void runByRow(Lhs& lhs, Rhs const& rhs )
  { for (int i = rhs.begin(); i< rhs.end(); ++i) { lhs.elt(i) = rhs.elt(i);}}
};


//---------------------NUMBER----------------------------------
//  point_ <- diagonal
template < typename Lhs, typename Rhs>
struct Copycat<  Lhs,  Rhs, Arrays::number_, Arrays::number_>
{
  inline static void runByCol(Lhs& lhs, Rhs const& rhs )
  { lhs.elt() = rhs.elt();}
  inline static void runByRow(Lhs& lhs, Rhs const& rhs )
  { lhs.elt() = rhs.elt();}
};

/** @ingroup hidden
 * utility class that select if the copy will be by row or by column
 **/
template < typename Lhs, typename Rhs, int TOrient_>
struct CopycatSelector;

/** specialization for column oriented arrrays */
template< typename Lhs, typename Rhs>
struct CopycatSelector< Lhs, Rhs, Arrays::by_col_>
{
  enum
  { tstructure_ = hidden::Traits<Lhs>::structure_
  , sstructure_ = hidden::Traits<Rhs>::structure_
  };
  inline static void run(Lhs& lhs, Rhs const& rhs )
  { Copycat<Lhs, Rhs, tstructure_, sstructure_>::runByCol(lhs, rhs );}
};

/** specialization for row oriented arrays */
template< typename Lhs, typename Rhs>
struct CopycatSelector< Lhs, Rhs, Arrays::by_row_>
{
  enum
  { tstructure_ = hidden::Traits<Lhs>::structure_
  , sstructure_ = hidden::Traits<Rhs>::structure_
  };
  inline static void run(Lhs& lhs, Rhs const& rhs )
  { Copycat<Lhs, Rhs, tstructure_, sstructure_>::runByRow(lhs, rhs );}
};

/** @ingroup hidden
 * utility class that select the resize method to call
 **/
template< typename Lhs, typename Rhs, int TStructure_>
struct resizeSelector;

/** 2D general case */
template< typename Lhs, typename Rhs, int TStructure_>
struct resizeSelector
{
  inline static void run(Lhs& lhs, ExprBase<Rhs> const& rhs )
  { lhs.resize(rhs.rows(), rhs.cols());}
};
/** specialization for the square_ case */
template< typename Lhs, typename Rhs>
struct resizeSelector<Lhs, Rhs, Arrays::square_>
{
  inline static void run(Lhs& lhs, ExprBase<Rhs> const& rhs )
  { lhs.resize(rhs.range());}
};
/** specialization for the diagonal_ case */
template< typename Lhs, typename Rhs>
struct resizeSelector<Lhs, Rhs, Arrays::diagonal_>
{
  inline static void run(Lhs& lhs, ExprBase<Rhs> const& rhs )
  { lhs.resize(rhs.range());}
};
/** specialization for the vector_ case */
template< typename Lhs, typename Rhs>
struct resizeSelector<Lhs, Rhs, Arrays::vector_>
{
  inline static void run(Lhs& lhs, ExprBase<Rhs> const& rhs )
  { lhs.resize(rhs.range());}
};
/** specialization for the point_ case */
template< typename Lhs, typename Rhs>
struct resizeSelector<Lhs, Rhs, Arrays::point_>
{
  inline static void run(Lhs& lhs, ExprBase<Rhs> const& rhs )
  { lhs.resize(rhs.range());}
};

/** specialization for the number_ case */
template< typename Lhs, typename Rhs>
struct resizeSelector<Lhs, Rhs, Arrays::number_>
{
  inline static void run(Lhs& lhs, ExprBase<Rhs> const& rhs )
  { /* nothing to do */;}
};

}  // namespace hidden


} // namespace STK

#endif /* STK_ASSIGNSELECTOR_H */
