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
// License for the specific language governing rights and limitations under
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

namespace Nektar
{
    typedef double NekDouble;

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
            std::map<FlagType, bool>::const_iterator x;
            return ((x = m_data.find(key)) != m_data.end() && x->second);
        }
    private:
        std::map<FlagType, bool> m_data;
    };

    /// An empty flag list.
    static FlagList NullFlagList;

} //end of namespace

#endif

/***
$Log: NektarUnivTypeDefs.hpp,v $
Revision 1.16  2008/06/10 06:00:37  bnelson
Updated documentation.

Revision 1.15  2008/05/29 21:32:11  pvos
Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions

Revision 1.14  2008/03/03 02:25:46  bnelson
Changed OneD, TwoD, and ThreeD to classes instead of enums to support type parameters in NekVector instead of unsigned int for the dimensions.

Revision 1.13  2007/08/10 03:37:21  jfrazier
Removed that cursed Mac formatting.

Revision 1.12  2007/08/06 05:35:45  ehan
Added enumerations for the 3 principal directions. These are to be used in GetD(Direction).

Revision 1.11  2007/07/22 23:03:24  bnelson
Backed out Nektar::ptr.

Revision 1.10  2007/07/20 00:39:54  bnelson
Replaced boost::shared_ptr with Nektar::ptr

Revision 1.9  2007/05/14 23:47:43  bnelson
Removed old SharedArray typedefs.  Added new Dimension enumeration.

Revision 1.8  2007/04/29 00:31:12  jfrazier
Continued conversion to multi_arrays.

Revision 1.7  2007/04/26 21:52:09  jfrazier
Converted to new multi_array implementation.

Revision 1.6  2007/04/10 14:00:44  sherwin
Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions

Revision 1.5  2007/03/31 15:37:51  bnelson
Removed boost::shared_array.hpp

Revision 1.4  2007/03/29 18:44:11  bnelson
Replaced boost::shared_array with SharedArray

Revision 1.3  2007/03/20 13:48:05  sherwin
Compiling version

Revision 1.2  2007/03/20 12:27:39  sherwin
Changed shared_ptr to shared_array

Revision 1.1  2007/03/20 11:56:25  sherwin
.

**/
