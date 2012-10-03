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
#include <LibUtilities/Foundations/FourierSingleModePoints.h>
#include <LibUtilities/Foundations/BLPoints.h>
#include <LibUtilities/Foundations/PolyEPoints.h>
#include <LibUtilities/Foundations/NodalTriElec.h>
#include <LibUtilities/Foundations/NodalTetElec.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/NodalTriEvenlySpaced.h>
#include <LibUtilities/Foundations/NodalTetEvenlySpaced.h>
#include <LibUtilities/Foundations/NodalPrismEvenlySpaced.h>
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
            const bool gaussInited1 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussLegendre), GaussPoints::Create);
            const bool gaussInited2 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMLegendre), GaussPoints::Create);
            const bool gaussInited3 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPLegendre), GaussPoints::Create);
            const bool gaussInited4 = PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoLegendre), GaussPoints::Create);
            const bool gaussInited5 = PointsManager().RegisterCreator(PointsKey(0, eGaussGaussChebyshev), GaussPoints::Create);
            const bool gaussInited6 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMChebyshev), GaussPoints::Create);
            const bool gaussInited7 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauPChebyshev), GaussPoints::Create);
            const bool gaussInited8 = PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoChebyshev), GaussPoints::Create);
            const bool gaussInited9 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta1), GaussPoints::Create);
            const bool gaussInited10 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha0Beta2), GaussPoints::Create);
            const bool gaussInited11 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha1Beta0), GaussPoints::Create);
            const bool gaussInited12 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauMAlpha2Beta0), GaussPoints::Create);
            
            const bool gaussInited13 = PointsManager().RegisterCreator(PointsKey(0, eGaussKronrodLegendre), GaussPoints::Create);
            const bool gaussInited14 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauKronrodMLegendre), GaussPoints::Create);
            const bool gaussInited15 = PointsManager().RegisterCreator(PointsKey(0, eGaussRadauKronrodMAlpha1Beta0), GaussPoints::Create);
            const bool gaussInited16 = PointsManager().RegisterCreator(PointsKey(0, eGaussLobattoKronrodLegendre), GaussPoints::Create);
        
            const bool fourierInited0 = PointsManager().RegisterCreator(PointsKey(0, eFourierEvenlySpaced), FourierPoints::Create);
            const bool fourierInitedSM0 = PointsManager().RegisterCreator(PointsKey(0, eFourierSingleModeSpaced), FourierSingleModePoints::Create);
            const bool BLInited0 = PointsManager().RegisterCreator(PointsKey(0, eBoundaryLayerPoints), BLPoints::Create);
            const bool polyeInited10 =  PointsManager().RegisterCreator(PointsKey(0, ePolyEvenlySpaced), PolyEPoints::Create);
            const bool NodalTriInited0 =  PointsManager().RegisterCreator(PointsKey(0, eNodalTriElec), NodalTriElec::Create);
            const bool NodalTriInited1 =  PointsManager().RegisterCreator(PointsKey(0, eNodalTriFekete), NodalTriFekete::Create);
            const bool nodalTetElecInited = PointsManager().RegisterCreator(PointsKey(0, eNodalTetElec), NodalTetElec::Create);
            const bool NodalTriEveInited = PointsManager().RegisterCreator(PointsKey(0, eNodalTriEvenlySpaced), NodalTriEvenlySpaced::Create);
            const bool NodalTetEveInited = PointsManager().RegisterCreator(PointsKey(0, eNodalTetEvenlySpaced), NodalTetEvenlySpaced::Create);
            const bool NodalPrismEveInited = PointsManager().RegisterCreator(PointsKey(0, eNodalPrismEvenlySpaced), NodalPrismEvenlySpaced::Create);

            const bool Basis_Inited = BasisManager().RegisterGlobalCreator(Basis::Create);
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


