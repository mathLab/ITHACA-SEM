#ifndef NEKTAR_LIBUTILITIES_BASICUTILS_STRINGUTILS_HPP
#define NEKTAR_LIBUTILITIES_BASICUTILS_STRINGUTILS_HPP

#include <string>
#include <cctype>
#include <algorithm>

namespace Nektar
{

inline void to_upper(std::string &s)
{
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) {
                       return std::tolower(c);
                   });
}

inline std::string to_upper_copy(const std::string &s)
{
    std::string cp = s;
    to_upper(cp);
    return cp;
}

}

#endif
