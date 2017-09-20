////////////////////////////////////////////////////////////////////////////////
//
//  File:  ParseUtils.hpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
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

#include <sstream>
#include <boost/spirit/include/qi.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;

namespace Nektar
{

/**
 * Helper class for boost::spirit. This pushes back values onto a
 * std::vector.
 */
template<typename T>
struct PushBackFunctor
{
    PushBackFunctor(std::vector<T> &in) : vec(in) {}

    void operator()(T num) const
    {
        vec.push_back(num);
    }
private:
    std::vector<T> &vec;
};

/**
 * Helper class for boost::spirit. This pushes back a range of values onto a
 * std::vector. Valid only for integer types.
 */
template<typename T>
struct SeqFunctor
{
    SeqFunctor(std::vector<T> &in) : vec(in) {}

    void operator()(fusion::vector<T, T> num) const
    {
        static_assert(std::is_integral<T>::value, "Integer type required.");
        for (int i = fusion::at_c<0>(num); i <= fusion::at_c<1>(num); ++i)
        {
            vec.push_back(i);
        }
    }
private:
    std::vector<T> &vec;
};

bool ParseUtils::GenerateSeqVector(
    const std::string &str, std::vector<unsigned int> &out)
{
    PushBackFunctor<unsigned int> f1(out);
    SeqFunctor<unsigned int> f2(out);

    return qi::phrase_parse(
        str.begin(),
        str.end(),
        ((qi::uint_ >> '-' >> qi::uint_)[f2] | qi::uint_[f1]) % ',',
        qi::ascii::space);
}

template <typename T>
bool ParseUtils::GenerateVector(const std::string &str, std::vector<T> &out)
{
    return qi::phrase_parse(
        str.begin(), str.end(), qi::auto_ % ',', qi::ascii::space, out);
}

template bool ParseUtils::GenerateVector<int>(
    const std::string &str, std::vector<int> &out);
template bool ParseUtils::GenerateVector<long>(
    const std::string &str, std::vector<long> &out);
template bool ParseUtils::GenerateVector<unsigned int>(
    const std::string &str, std::vector<unsigned int> &out);
template bool ParseUtils::GenerateVector<double>(
    const std::string &str, std::vector<double> &out);
template bool ParseUtils::GenerateVector<float>(
    const std::string &str, std::vector<float> &out);

template <>
bool ParseUtils::GenerateVector(const std::string &str,
                                std::vector<std::string> &out)
{
    return qi::phrase_parse(
        str.begin(), str.end(), *~qi::char_(",") % ',', qi::ascii::space, out);
}

}
