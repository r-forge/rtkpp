/*
 * STK_Kernel_Util.h
 *
 *  Created on: Nov 22, 2017
 *      Author: iovleff
 */

#ifndef STK_KERNEL_UTIL_H
#define STK_KERNEL_UTIL_H

namespace STK
{

namespace hidden
{
/** @ingroup hidden, STKernel
 *  @brief The traits struct IteratorTraits must be specialized for any iterator
 *  derived from the base class IteratorBase.
 *  We use the type names defined by the STL for the iterator_traits class.
 *
 *  For example:
 *  @code
 *  template<typename Type>
 *  struct IteratorTraits
 *  {
 *    /// One of the iterator_tags types
 *    typedef std::random_access_iterator_tag  iterator_category;
 *    /// The type "pointed to" by the iterator.
 *    typedef Type        value_type;
 *    /// Distance between iterators is represented as this type.
 *    typedef int  difference_type;
 *    /// This type represents a pointer-to-value_type.
 *    typedef Type*   pointer;
 *    /// This type represents a reference-to-value_type.
 *    typedef Type& reference;
 *  };
 *  @endcode
 */
template <typename Derived> struct IteratorTraits;

} // namespace hidden

} // namespace STK


#endif /* STK_KERNEL_UTIL_H_ */
