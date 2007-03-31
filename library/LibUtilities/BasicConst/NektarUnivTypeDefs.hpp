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

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <vector>

namespace Nektar
{

    typedef double NekDouble;
    typedef SharedArray<NekDouble> NekDoubleSharedArray;
    typedef std::vector< NekDoubleSharedArray > NekDoubleArrayVector;
    typedef std::vector< NekDoubleSharedArray >::iterator NekDoubleArrayVectorIter;
 
    typedef SharedArray<int>  NekIntSharedArray;
    typedef std::vector< NekIntSharedArray > NekIntArrayVector;
    typedef std::vector< NekIntSharedArray >::iterator NekIntArrayVectorIter;

} //end of namespace 

#endif

/***
$Log: NektarUnivTypeDefs.hpp,v $
Revision 1.4  2007/03/29 18:44:11  bnelson
Replaced boost::shared_array with SharedArray

Revision 1.3  2007/03/20 13:48:05  sherwin
Compiling version

Revision 1.2  2007/03/20 12:27:39  sherwin
Changed shared_ptr to shared_array

Revision 1.1  2007/03/20 11:56:25  sherwin
.

**/
