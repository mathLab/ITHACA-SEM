#pragma once

#include <type_traits>
#include <boost/align/aligned_allocator.hpp>
#include "traits.hpp"

namespace tinysimd
{

// should be enabled only for vector types
template <typename T, typename = typename std::enable_if<is_vector<T>::value>::type>
using allocator = boost::alignment::aligned_allocator<T, T::alignment>;

} // namespace tinysimd

