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

#include <boost/multi_array.hpp>
#include <LibUtilities/BasicUtils/SharedPtr.hpp>
#include <vector>

namespace Nektar
{
    typedef double NekDouble;

    enum Dimension
    {
        OneD = 1,
        TwoD = 2,
        ThreeD = 3
    };

} //end of namespace 

#endif

/***
$Log: NektarUnivTypeDefs.hpp,v $
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
