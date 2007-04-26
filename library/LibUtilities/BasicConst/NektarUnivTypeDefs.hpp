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
#include <boost/shared_ptr.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <vector>

namespace Nektar
{
    typedef double NekDouble;

    // Multi_array and multi_array_ref type defines, and their
    // shared_ptr versions.  The ref version is used to create
    // a multi_array that uses memory not owned by the multi_array.

    template<typename T>
    struct Nek1DArray
    {
        typedef boost::multi_array<T, 1> type;
    };

    template<typename T>
    struct Nek1DSharedArray
    {
        typedef boost::shared_ptr<typename Nek1DArray<T>::type> type;
    };

    template<typename T>
    struct Nek1DConstSharedArray
    {
        typedef boost::shared_ptr<const typename Nek1DArray<T>::type> type;
    };

    typedef Nek1DArray<NekDouble>::type                 NekDouble1DArray;
    typedef Nek1DSharedArray<NekDouble>::type           NekDouble1DSharedArray;
    typedef Nek1DConstSharedArray<NekDouble>::type      ConstNekDouble1DSharedArray;
    typedef Nek1DArray<int>::type                       Int1DArray;
    typedef Nek1DSharedArray<int>::type                 Int1DSharedArray;
    typedef Nek1DConstSharedArray<int>::type            ConstInt1DSharedArray;

    template<typename T>
    struct Nek2DArray
    {
        typedef boost::multi_array<T, 2> type;
    };

    template<typename T>
    struct Nek2DSharedArray
    {
        typedef boost::shared_ptr<typename Nek2DArray<T>::type> type;
    };

    template<typename T>
    struct Nek2DConstSharedArray
    {
        typedef boost::shared_ptr<const typename Nek2DArray<T>::type> type;
    };

    typedef Nek2DArray<NekDouble>::type                 NekDouble2DArray;
    typedef Nek2DSharedArray<NekDouble>::type           NekDouble2DSharedArray;
    typedef Nek2DConstSharedArray<NekDouble>::type      ConstNekDouble2DSharedArray;
    typedef Nek2DArray<int>::type                       Int2DArray;
    typedef Nek2DSharedArray<int>::type                 Int2DSharedArray;
    typedef Nek2DConstSharedArray<int>::type            ConstInt2DSharedArray;

    template<typename T>
    struct Nek3DArray
    {
        typedef boost::multi_array<T, 3> type;
    };

    template<typename T>
    struct Nek3DSharedArray
    {
        typedef boost::shared_ptr<typename Nek3DArray<T>::type> type;
    };

    template<typename T>
    struct Nek3DConstSharedArray
    {
        typedef boost::shared_ptr<const typename Nek3DArray<T>::type> type;
    };

    typedef Nek3DArray<NekDouble>::type                 NekDouble3DArray;
    typedef Nek3DSharedArray<NekDouble>::type           NekDouble3DSharedArray;
    typedef Nek3DConstSharedArray<NekDouble>::type      ConstNekDouble3DSharedArray;
    typedef Nek3DArray<int>::type                       Int3DArray;
    typedef Nek3DSharedArray<int>::type                 Int3DSharedArray;
    typedef Nek3DConstSharedArray<int>::type            ConstInt3DSharedArray;

    // Old typedefs based on SharedArray.  Easiest way to force the issue is to
    // comment them out and then recompile.

    //typedef std::vector< NekDouble1DSharedArray > NekDoubleArrayVector;
    //typedef std::vector< NekDouble1DSharedArray >::iterator NekDoubleArrayVectorIter;
 
    //typedef SharedArray<int>  NekIntSharedArray;
    //typedef std::vector< NekIntSharedArray > NekIntArrayVector;
    //typedef std::vector< NekIntSharedArray >::iterator NekIntArrayVectorIter;

    //typedef SharedArray<const int>  ConstNekIntSharedArray;

} //end of namespace 

#endif

/***
$Log: NektarUnivTypeDefs.hpp,v $
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
