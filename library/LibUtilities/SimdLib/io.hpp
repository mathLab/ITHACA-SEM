#pragma once

#include "traits.hpp"
#include "tinysimd.hpp"
#include <ostream>

namespace tinysimd
{


template <class T, typename = typename std::enable_if
        <tinysimd::is_vector<T>::value>::type
    >
std::ostream& operator<<(std::ostream& os, const T& avec)
{
    alignas(T::alignment) typename T::scalarArray tmp;
    avec.store(tmp);
    for (unsigned short i = 0; i < T::width; ++i)
    {
        os << tmp[i] << '\t';
    }
    return os;
}


} // namespace tinysimd