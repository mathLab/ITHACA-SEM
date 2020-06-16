#pragma once

#include <type_traits>
#include <boost/align/aligned_allocator.hpp>

namespace tinysimd
{

// should be enabled only for vector types
template<typename T>
using allocator = boost::alignment::aligned_allocator<T, T::alignment>;

} // namespace tinysimd

