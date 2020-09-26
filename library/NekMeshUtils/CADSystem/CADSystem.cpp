////////////////////////////////////////////////////////////////////////////////
//
//  File: CADSystem.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: cad object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/CADSystem/CADSystem.h>

using namespace std;

namespace Nektar
{
namespace NekMeshUtils
{

EngineFactory &GetEngineFactory()
{
    static EngineFactory instance;
    return instance;
}

CADVertFactory &GetCADVertFactory()
{
    static CADVertFactory instance;
    return instance;
}

CADCurveFactory &GetCADCurveFactory()
{
    static CADCurveFactory instance;
    return instance;
}

CADSurfFactory &GetCADSurfFactory()
{
    static CADSurfFactory instance;
    return instance;
}

Array<OneD, NekDouble> CADSystem::GetPeriodicTranslationVector(int first,
                                                               int second)
{
    ASSERTL0(GetNumSurf() == 1, "wont work for multi surfaces yet");

    CADCurveSharedPtr c1 = GetCurve(first);
    CADCurveSharedPtr c2 = GetCurve(second);

    NekDouble tst = c1->GetTotLength() - c2->GetTotLength();
    ASSERTL0(fabs(tst) < 1e-6, "periodic curves not same length");

    vector<CADVertSharedPtr> v1 = c1->GetVertex();
    Array<OneD, NekDouble> p1 = v1[0]->GetLoc();

    Array<OneD, NekDouble> p2;
    vector<CADVertSharedPtr> v2 = c2->GetVertex();
    if (c1->GetOrienationWRT(1) == c2->GetOrienationWRT(1))
    {
        p2 = v2[1]->GetLoc();
    }
    else
    {
        p2 = v2[0]->GetLoc();
    }

    Array<OneD, NekDouble> ret(3);
    ret[0] = p2[0] - p1[0];
    ret[1] = p2[1] - p1[1];
    ret[2] = p2[2] - p1[2];

    return ret;
}

}
}
