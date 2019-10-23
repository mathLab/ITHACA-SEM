////////////////////////////////////////////////////////////////////////////////
//
//  File:  ParseUtils.cpp
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

#include <sstream>
#include <boost/core/ignore_unused.hpp>
#include <boost/spirit/include/qi_core.hpp>
#include <boost/spirit/include/qi_auto.hpp>
#include <LibUtilities/BasicUtils/ParseUtils.h>

namespace qi = boost::spirit::qi;
namespace fusion = boost::fusion;

namespace Nektar
{

/**
 * @brief Helper functors for holding a vector of numbers to be parsed by
 * boost::spirit.
 *
 * @see ParseUtils::GenerateSeqVector
 */
template<typename T>
struct PushBackFunctor
{
    PushBackFunctor(std::vector<T> &in) : m_vec(in) {}

    /**
     * @brief Pushes back values onto #m_vec as given by @p num.
     */
    void operator()(T num) const
    {
        m_vec.push_back(num);
    }

    /**
     * @brief Pushes back values onto #m_vec between the range supplied by @p
     * num. Valid for only integer types.
     */
    void operator()(fusion::vector<T, T> num) const
    {
        static_assert(std::is_integral<T>::value, "Integer type required.");
        for (T i = fusion::at_c<0>(num); i <= fusion::at_c<1>(num); ++i)
        {
            m_vec.push_back(i);
        }
    }
private:
    // Do not allow assignment
    PushBackFunctor& operator=(const PushBackFunctor& src)
    {
        boost::ignore_unused(src);
        return *this;
    }

    /// Storage vector that will hold parsed variables from boost::spirit.
    std::vector<T> &m_vec;
};

/**
 * @brief Takes a comma-separated compressed string and converts it to entries
 * in a vector.
 *
 * This routine is the inverse of ParseUtils::GenerateSeqString. For example,
 *
 *     std::string input = "1-4,6-8,5,2,3";
 *     std::vector<unsigned int> output;
 *     ParseUtils::GenerateSeqString(input, output);
 *
 * produces an `output` vector with the entries `{1,2,3,4,6,7,8,5,2,3}`.
 *
 * @param str  Input CSV string of unsigned integers.
 * @param out  Output vector.
 *
 * @see ParseUtils::GenerateSeqString
 */
bool ParseUtils::GenerateSeqVector(
    const std::string &str, std::vector<unsigned int> &out)
{
    PushBackFunctor<unsigned int> f1(out), f2(out);

    auto it = str.begin();
    bool success = qi::phrase_parse(
        it,
        str.end(),
        ((qi::uint_ >> '-' >> qi::uint_)[f2] | qi::uint_[f1]) % ',',
        qi::ascii::space);

    return success && it == str.end();
}

/**
 * @brief Takes a comma-separated string and converts it to entries in a vector.
 *
 * This routine splits up a comma-separated string and returns a vector with the
 * entries. Template specialisations should be defined in this file (and not in
 * the header file) as the use of boost::spirit::qi makes compilation times
 * quite slow.
 *
 * @param str  Input CSV string.
 * @param out  Output vector.
 */
template <typename T>
bool ParseUtils::GenerateVector(const std::string &str, std::vector<T> &out)
{
    auto it = str.begin();
    bool success = qi::phrase_parse(
        it, str.end(), qi::auto_ % ',', qi::ascii::space, out);
    return success && it == str.end();
}

template LIB_UTILITIES_EXPORT bool ParseUtils::GenerateVector<int>(
    const std::string &str, std::vector<int> &out);
template LIB_UTILITIES_EXPORT bool ParseUtils::GenerateVector<long>(
    const std::string &str, std::vector<long> &out);
template LIB_UTILITIES_EXPORT bool ParseUtils::GenerateVector<unsigned int>(
    const std::string &str, std::vector<unsigned int> &out);
template LIB_UTILITIES_EXPORT bool ParseUtils::GenerateVector<double>(
    const std::string &str, std::vector<double> &out);
template LIB_UTILITIES_EXPORT bool ParseUtils::GenerateVector<float>(
    const std::string &str, std::vector<float> &out);

/**
 * @brief Specialised version of ParseUtils::GenerateVector for std::string.
 *
 * This routine specialises for the std::string data type as this type is not
 * supported by boost::spirit::qi::auto_.
 */
template <> LIB_UTILITIES_EXPORT
bool ParseUtils::GenerateVector(const std::string &str,
                                std::vector<std::string> &out)
{
    auto it = str.begin();
    bool success = qi::phrase_parse(
        it, str.end(), +~qi::char_(",") % ',', qi::ascii::space, out);
    return success && it == str.end();
}

}
