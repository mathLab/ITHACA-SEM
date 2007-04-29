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

#include <iostream>
#include <loki/Singleton.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/FourierPoints.h>
#include <LibUtilities/Foundations/PolyEPoints.h>
#include <LibUtilities/Foundations/NodalTriElec.h>
#include <LibUtilities/Foundations/NodalTetElec.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>

namespace Nektar
{
    namespace LibUtilities 
    {
	// Register all points and basis creators.
	namespace
        {
            const bool gaussInited1 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussLegendre), GaussPoints<NekDouble>::Create);
            const bool gaussInited2 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMLegendre), GaussPoints<NekDouble>::Create);
            const bool gaussInited3 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPLegendre), GaussPoints<NekDouble>::Create);
            const bool gaussInited4 = PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoLegendre), GaussPoints<NekDouble>::Create);
            const bool gaussInited5 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussChebyshev), GaussPoints<NekDouble>::Create);
            const bool gaussInited6 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMChebyshev), GaussPoints<NekDouble>::Create);
            const bool gaussInited7 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPChebyshev), GaussPoints<NekDouble>::Create);
            const bool gaussInited8 = PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoChebyshev), GaussPoints<NekDouble>::Create);
            const bool gaussInited9 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1), GaussPoints<NekDouble>::Create);
            const bool gaussInited10 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2), GaussPoints<NekDouble>::Create);
	    
            const bool fourierInited0 = PointsManager().RegisterCreator(PointsKey(0, eFourierEvenlySpaced), FourierPoints<NekDouble>::Create);
            const bool polyeInited10 =  PointsManager().RegisterCreator(PointsKey(0, ePolyEvenlySpaced), PolyEPoints<NekDouble>::Create);
            const bool NodalTriInited0 =  PointsManager().RegisterCreator(PointsKey(0, eNodalTriElec), NodalTriElec<NekDouble>::Create);
            const bool NodalTriInited1 =  PointsManager().RegisterCreator(PointsKey(0, eNodalTriFekete), NodalTriFekete<NekDouble>::Create);
            const bool nodalTetElecInited = PointsManager().RegisterCreator(PointsKey(0, eNodalTetElec), NodalTetElec<NekDouble>::Create);

	        const bool Ortho_A_Inited = BasisManager().RegisterCreator(BasisKey(eOrtho_A, 0, NullPointsKey), Basis::Create);
            const bool Ortho_B_Inited = BasisManager().RegisterCreator(BasisKey(eOrtho_B, 0, NullPointsKey), Basis::Create);
            const bool Ortho_C_Inited = BasisManager().RegisterCreator(BasisKey(eOrtho_C, 0, NullPointsKey), Basis::Create);
            const bool Modified_A_Inited = BasisManager().RegisterCreator(BasisKey(eModified_A, 0, NullPointsKey), Basis::Create);
            const bool Modified_B_Inited = BasisManager().RegisterCreator(BasisKey(eModified_B, 0, NullPointsKey), Basis::Create);
            const bool Modified_C_Inited = BasisManager().RegisterCreator(BasisKey(eModified_C, 0, NullPointsKey), Basis::Create);
            const bool Fourier_Inited = BasisManager().RegisterCreator(BasisKey(eFourier, 0, NullPointsKey), Basis::Create);
            const bool GLL_Lagrange_Inited = BasisManager().RegisterCreator(BasisKey(eGLL_Lagrange, 0, NullPointsKey), Basis::Create);
            const bool Legendre_Inited = BasisManager().RegisterCreator(BasisKey(eLegendre, 0, NullPointsKey), Basis::Create);
            const bool Chebyshev_Inited = BasisManager().RegisterCreator(BasisKey(eChebyshev, 0, NullPointsKey), Basis::Create);
        };

        PointsManagerT &PointsManager(void)
        {
            return Loki::SingletonHolder<PointsManagerT>::Instance();
        }

        BasisManagerT &BasisManager(void)
        {
            return Loki::SingletonHolder<BasisManagerT>::Instance();
        }

    } // end of namespace LibUtilities
} // end of namespace Nektar


/**
$Log: ManagerAccess.cpp,v $
Revision 1.10  2007/04/10 14:00:45  sherwin
Update to include SharedArray in all 2D element (including Nodal tris). Have also remvoed all new and double from 2D shapes in StdRegions

Revision 1.9  2007/04/10 02:40:53  bnelson
Updated loki capitalization.

Revision 1.8  2007/03/19 12:46:16  kirby
*** empty log message ***

Revision 1.7  2007/03/13 21:31:32  kirby
 small update to the numbering.

Revision 1.6  2007/03/13 16:59:04  kirby
 added GaussLobattoChebyshev -- we had forgotten it

Revision 1.5  2007/02/26 15:52:31  sherwin
Working version for Fourier points calling from StdRegions. Fourier interpolations not working yet

Revision 1.4  2007/02/06 17:12:31  jfrazier
Fixed a problem with global initialization in libraries.

Revision 1.3  2007/02/01 23:28:42  jfrazier
Basis is not working, but not fully tested.

Revision 1.2  2007/01/20 21:45:59  kirby
*** empty log message ***

Revision 1.1  2007/01/19 18:02:26  jfrazier
Initial checkin.

**/



