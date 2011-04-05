///////////////////////////////////////////////////////////////////////////////
//
// File Points1D.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
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
// Description: C functions to provide access to managers. 
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILIITIES_FOUNDATIONS_MANAGER_ACCESS_H
#define NEKTAR_LIB_UTILIITIES_FOUNDATIONS_MANAGER_ACCESS_H

#include <LibUtilities/BasicUtils/NekManager.hpp>
#include <LibUtilities/Foundations/Points.h>
#include <LibUtilities/Foundations/Basis.h>

namespace Nektar
{
    namespace LibUtilities 
    {
        typedef NekManager<PointsKey, Points<NekDouble>, PointsKey::opLess> PointsManagerT;
        LIB_UTILITIES_EXPORT PointsManagerT &PointsManager(void);

        typedef NekManager<BasisKey, Basis, BasisKey::opLess> BasisManagerT;
        LIB_UTILITIES_EXPORT BasisManagerT &BasisManager(void);

    } // end of namespace LibUtilities
} // end of namespace Nektar
#endif //NEKTAR_LIB_UTILIITIES_FOUNDATIONS_MANAGER_ACCESS_H

/**
$Log: ManagerAccess.h,v $
Revision 1.10  2008/07/12 11:37:53  pvos
Added time integration scheme manager

Revision 1.9  2007/04/29 03:09:47  jfrazier
More conversion to multi_arrays.

Revision 1.8  2007/02/06 17:12:31  jfrazier
Fixed a problem with global initialization in libraries.

Revision 1.7  2007/02/01 23:28:42  jfrazier
Basis is working, but not fully tested.

Revision 1.6  2007/01/25 21:31:46  jfrazier
Format change.

Revision 1.5  2007/01/20 21:52:34  sherwin
Remove Basis template class definitino

Revision 1.4  2007/01/20 21:45:59  kirby
*** empty log message ***

Revision 1.3  2007/01/20 21:33:58  sherwin
Added method

Revision 1.2  2007/01/19 21:59:27  sherwin
Some SJS mods - still does not compile yet

Revision 1.1  2007/01/19 18:02:26  jfrazier
Initial checkin.

**/

