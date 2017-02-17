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

#include <LibUtilities/Foundations/BLPoints.h>
#include <LibUtilities/Foundations/Basis.h>
#include <LibUtilities/Foundations/Foundations.hpp>
#include <LibUtilities/Foundations/FourierPoints.h>
#include <LibUtilities/Foundations/FourierSingleModePoints.h>
#include <LibUtilities/Foundations/GaussPoints.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/Foundations/NodalHexElec.h>
#include <LibUtilities/Foundations/NodalPrismElec.h>
#include <LibUtilities/Foundations/NodalPrismEvenlySpaced.h>
#include <LibUtilities/Foundations/NodalPrismSPI.h>
#include <LibUtilities/Foundations/NodalQuadElec.h>
#include <LibUtilities/Foundations/NodalTetElec.h>
#include <LibUtilities/Foundations/NodalTetEvenlySpaced.h>
#include <LibUtilities/Foundations/NodalTetSPI.h>
#include <LibUtilities/Foundations/NodalTriElec.h>
#include <LibUtilities/Foundations/NodalTriEvenlySpaced.h>
#include <LibUtilities/Foundations/NodalTriFekete.h>
#include <LibUtilities/Foundations/NodalTriSPI.h>
#include <LibUtilities/Foundations/PolyEPoints.h>
#include <loki/Singleton.h>

namespace Nektar
{
namespace LibUtilities
{
// Register all points and basis creators.
namespace
{
// this is a trick class to silence the warnings from the RegisterCreators
class setup
{
public:
    static const bool gaussInited1, gaussInited2, gaussInited3, gaussInited4,
        gaussInited5, gaussInited6, gaussInited7, gaussInited8, gaussInited9,
        gaussInited10, gaussInited11, gaussInited12, gaussInited13,
        gaussInited14, gaussInited15, gaussInited16;
    static const bool fourierInited1, fourierInited2;
    static const bool BLInited1, BLInited2;
    static const bool polyeInited1;
    static const bool NodalTriInited1, NodalTriInited2, NodalTriInited3,
        NodalTriInited4;
    static const bool NodalPrismInited1, NodalPrismInited2, NodalPrismInited3;
    static const bool NodalQuadInited1;
    static const bool NodalTetInited1, NodalTetInited2, NodalTetInited3;
    static const bool NodalHexInited1;
    static const bool Basis_Inited;
};

const bool setup::gaussInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussGaussLegendre), GaussPoints::Create);
const bool setup::gaussInited2 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauMLegendre), GaussPoints::Create);
const bool setup::gaussInited3 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauPLegendre), GaussPoints::Create);
const bool setup::gaussInited4 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussLobattoLegendre), GaussPoints::Create);
const bool setup::gaussInited5 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussGaussChebyshev), GaussPoints::Create);
const bool setup::gaussInited6 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauMChebyshev), GaussPoints::Create);
const bool setup::gaussInited7 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauPChebyshev), GaussPoints::Create);
const bool setup::gaussInited8 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussLobattoChebyshev), GaussPoints::Create);
const bool setup::gaussInited9 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauMAlpha0Beta1), GaussPoints::Create);
const bool setup::gaussInited10 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauMAlpha0Beta2), GaussPoints::Create);
const bool setup::gaussInited11 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauMAlpha1Beta0), GaussPoints::Create);
const bool setup::gaussInited12 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauMAlpha2Beta0), GaussPoints::Create);
const bool setup::gaussInited13 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussKronrodLegendre), GaussPoints::Create);
const bool setup::gaussInited14 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauKronrodMLegendre), GaussPoints::Create);
const bool setup::gaussInited15 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussRadauKronrodMAlpha1Beta0), GaussPoints::Create);
const bool setup::gaussInited16 = PointsManager().RegisterCreator(
    PointsKey(0, eGaussLobattoKronrodLegendre), GaussPoints::Create);

const bool setup::fourierInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eFourierEvenlySpaced), FourierPoints::Create);
const bool setup::fourierInited2 = PointsManager().RegisterCreator(
    PointsKey(0, eFourierSingleModeSpaced), FourierSingleModePoints::Create);

const bool setup::BLInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eBoundaryLayerPoints), BLPoints::Create);
const bool setup::BLInited2 = PointsManager().RegisterCreator(
    PointsKey(0, eBoundaryLayerPointsRev), BLPoints::Create);

const bool setup::polyeInited1 = PointsManager().RegisterCreator(
    PointsKey(0, ePolyEvenlySpaced), PolyEPoints::Create);

const bool setup::NodalTriInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalTriElec), NodalTriElec::Create);
const bool setup::NodalTriInited2 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalTriFekete), NodalTriFekete::Create);
const bool setup::NodalTriInited3 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalTriSPI), NodalTriSPI::Create);
const bool setup::NodalTriInited4 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalTriEvenlySpaced), NodalTriEvenlySpaced::Create);

const bool setup::NodalQuadInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalQuadElec), NodalQuadElec::Create);

const bool setup::NodalTetInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalTetElec), NodalTetElec::Create);
const bool setup::NodalTetInited2 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalTetSPI), NodalTetSPI::Create);
const bool setup::NodalTetInited3 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalTetEvenlySpaced), NodalTetEvenlySpaced::Create);

const bool setup::NodalPrismInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalPrismEvenlySpaced), NodalPrismEvenlySpaced::Create);
const bool setup::NodalPrismInited2 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalPrismElec), NodalPrismElec::Create);
const bool setup::NodalPrismInited3 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalPrismSPI), NodalPrismSPI::Create);

const bool setup::NodalHexInited1 = PointsManager().RegisterCreator(
    PointsKey(0, eNodalHexElec), NodalHexElec::Create);

const bool setup::Basis_Inited =
    BasisManager().RegisterGlobalCreator(Basis::Create);
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
