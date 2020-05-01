///////////////////////////////////////////////////////////////////////////////
//
// File: InterpBasis.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Demo for testing functionality of StdProject
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <LibUtilities/BasicUtils/Timer.h>

#include "StdDemoSupport.hpp"

int main(int argc, char *argv[])
{
    DemoSupport demo;
    demo.ParseArguments(argc, argv);
    StdExpansion *E = demo.CreateStdExpansion();

    if (E == nullptr)
    {
        return 1;
    }

    int nCoeffs = E->GetNcoeffs(), nPts = E->GetTotPoints();
    int nTot = nCoeffs * nPts, dimension = E->GetShapeDimension();

    Array<OneD, Array<OneD, NekDouble>> coords = demo.GetCoords(E);
    Array<OneD, NekDouble> sol(nTot), phys(nTot), tmpIn(dimension);

    // For each mode, we follow two approaches:
    //
    // 1) Evaluate the basis function at each quadrature point using the
    //    StdExpansion::PhysEvaluateBasis function.
    // 2) Evaluate the basis function at all quadrature points using FillMode.
    //
    // These are then compared to ensure they give the same result.
    for (int k = 0; k < nCoeffs; ++k)
    {
        // Evaluate each mode at the quadrature points.
        for (int i = 0; i < nPts; ++i)
        {
            for (int d = 0; d < dimension; ++d)
            {
                tmpIn[d] = coords[d][i];
            }

            phys[k * nPts + i] = E->PhysEvaluateBasis(tmpIn, k);
        }

        // Fill the 'solution' field with each of the modes using FillMode.
        Array<OneD, NekDouble> tmp = sol + k * nPts;
        E->FillMode(k, tmp);
    }

    cout << "L infinity error : " << scientific << E->Linf(phys, sol) << endl;
    cout << "L 2 error        : " << scientific << E->L2(phys, sol) << endl;

    return 0;
}
