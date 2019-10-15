///////////////////////////////////////////////////////////////////////////////
//
// File SessionReader.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Hash utilities for C++11 STL maps to support combined hashes
//              and enumerations.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef LIBUTILITIES_BASICUTILS_HASHUTILS
#define LIBUTILITIES_BASICUTILS_HASHUTILS

#include <boost/core/ignore_unused.hpp>

#include <functional>

namespace Nektar
{

inline void hash_combine(std::size_t& seed)
{
    boost::ignore_unused(seed);
}

template <typename T, typename... Args>
inline void hash_combine(std::size_t& seed, const T& v, Args... args)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    hash_combine(seed, args...);
}

template <typename T, typename... Args>
inline std::size_t hash_combine(const T& v, Args... args)
{
    boost::ignore_unused(v);
    std::size_t seed = 0;
    hash_combine(seed, args...);
    return seed;
}

template<typename Iter>
std::size_t hash_range(Iter first, Iter last)
{
    std::size_t seed = 0;
    for(; first != last; ++first)
    {
        hash_combine(seed, *first);
    }
    return seed;
}

template<typename Iter>
void hash_range(std::size_t &seed, Iter first, Iter last)
{
    hash_combine(seed, hash_range(first, last));
}

struct EnumHash
{
    template <typename T>
    std::size_t operator()(T t) const
    {
        return static_cast<std::size_t>(t);
    }
};

struct PairHash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const
    {
        std::size_t seed = 0;
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        hash_combine(seed, h1, h2);
        return seed;
    }
};

}
#endif
