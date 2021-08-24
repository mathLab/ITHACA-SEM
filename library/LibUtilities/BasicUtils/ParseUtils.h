////////////////////////////////////////////////////////////////////////////////
//
//  File:  ParseUtils.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:  This file contains various parsing utilities, primarily used
//                by SpatialDomains to process input files.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBUTILITIES_PARSEUTILS_H
#define NEKTAR_LIBUTILITIES_PARSEUTILS_H

#include <sstream>
#include <vector>
#include <LibUtilities/LibUtilitiesDeclspec.h>

namespace Nektar
{

class ParseUtils
{
public:
    LIB_UTILITIES_EXPORT static bool GenerateSeqVector(
        const std::string &str, std::vector<unsigned int> &out);

    template <typename T>
    static bool GenerateVector(const std::string &str, std::vector<T> &out);

    /**
     * @brief Generate a compressed comma-separated string representation of a
     * vector of unsigned integers.
     *
     * This utility routine takes entries of @p v and returns a string
     * sequence. For example,
     *
     *     std::vector<unsigned int> vec = {1,2,3,4,6,7,8,5,2,3};
     *     std::string output = ParseUtils::GenerateVector(vec);
     *
     * will produce an `output` string containing `1-4,6-8,5,2,3`.
     *
     * @param v  Vector of unsigned integers.
     * @return   Compressed comma separated string.
     */
    template <typename T>
    static std::string GenerateSeqString(const std::vector<T> &v)
    {
        static_assert(std::is_integral<T>::value && std::is_unsigned<T>::value,
                      "Unsigned integer type required.");

        if (v.size() == 0)
        {
            return "";
        }

        std::ostringstream ss;
        auto first = v[0], last = v[0];

        ss << v[0];

        for (auto &i : v)
        {
            if (i != last + 1 && i != last)
            {
                if (last != first)
                {
                    ss << '-' << last;
                }
                ss << ',' << i;
                first = i;
            }
            last = i;
        }

        if (last != first)
        {
            ss << '-' << last;
        }

        return ss.str();
    }
};

}

#endif //NEKTAR_LIBUTILITIES_PARSEUTILS_HPP
