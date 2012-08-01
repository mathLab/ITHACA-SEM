///////////////////////////////////////////////////////////////////////////////
//
// File NektarUnivConsts.hpp
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
// Description: Universal constants in the Nektar Library
//
///////////////////////////////////////////////////////////////////////////////

#ifndef  NEKTARUNIVCONSTS_HPP
#define  NEKTARUNIVCONSTS_HPP

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>

namespace Nektar
{
    namespace NekConstants
    {
        static const NekDouble kNekUnsetDouble = -9999;
        static const NekDouble kVertexTheSameDouble  = 1.0e-8;
        static const NekDouble kGeomFactorsTol = 1.0e-8;
        static const NekDouble kNekZeroTol = 1.0e-12;
        static const NekDouble kGeomRightAngleTol = 1e-14;
        static const NekDouble kNekSqrtTol = 1.0e-16;
        static const NekDouble kNekIterativeTol = 1e-09;
    }
} //end of namespace

#endif

/***
$Log: NektarUnivConsts.hpp,v $
Revision 1.8  2008/09/09 13:59:51  sherwin
Added NekZeroTol

Revision 1.7  2007/11/29 17:00:13  sherwin
Update to do with MultiRegions stuff

Revision 1.6  2007/05/14 23:25:15  bnelson
Removed unneeded code.

Revision 1.5  2007/05/14 23:24:40  bnelson
Removed unneeded code.

Revision 1.4  2007/04/26 21:51:54  jfrazier
Converted to new multi_array implementation.

Revision 1.3  2007/03/31 00:39:51  bnelson
*** empty log message ***

Revision 1.2  2007/03/21 20:56:42  sherwin
Update to change BasisSharedVector to boost::shared_array<BasisSharedPtr> and removed tthe Vector definitions in GetCoords and PhysDeriv

Revision 1.1  2007/03/20 11:56:25  sherwin
.

**/
