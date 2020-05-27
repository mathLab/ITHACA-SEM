///////////////////////////////////////////////////////////////////////////////
//
// File NektarUnivTypeDefs.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Scientific Computing and Imaging Institute,
// University of Utah (USA) and Department of Aeronautics, Imperial
// College London (UK).
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
// Description: Universal type defines  in the Nektar Library
//
///////////////////////////////////////////////////////////////////////////////

#ifndef  NEKTARUNIVTYPEDEF_HPP
#define  NEKTARUNIVTYPEDEF_HPP

#include <map>
#include <cstdint>

namespace Nektar
{
    typedef double NekDouble;

    typedef std::int32_t  NekInt;
    typedef std::int32_t  NekInt32;
    typedef std::int64_t  NekInt64;
    typedef std::uint32_t NekUInt;
    typedef std::uint32_t NekUInt32;
    typedef std::uint64_t NekUInt64;

    struct OneD
    {
        static const unsigned int Value = 1;
    };
    struct TwoD
    {
        static const unsigned int Value = 2;
    };
    struct ThreeD
    {
        static const unsigned int Value = 3;
    };

    struct FourD
    {
        static const unsigned int Value = 4;
    };

    enum Direction
    {
        xDir = 0,
        yDir = 1,
        zDir = 2
    };

    enum OutputFormat
    {
        eTecplot,
        eTecplotZones,
        eTecplotSingleBlock,
        eGmsh,
        eGnuplot
    };

    /// Enumeration of flags for passing a list of options.
    enum FlagType
    {
        eUseGlobal
    };

    /// String map for FlagType enumeration.
    const char* const FlagTypeMap[] = {
        "UseGlobal"
    };

    /// Defines a list of flags.
    class FlagList
    {
    public:
        void set(const FlagType &key, bool value)
        {
            m_data[key] = value;
        }
        bool isSet(const FlagType &key) const
        {
            auto x = m_data.find(key);
            return x != m_data.end() && x->second;
        }
    private:
        std::map<FlagType, bool> m_data;
    };

} //end of namespace

#endif
