///////////////////////////////////////////////////////////////////////////////
//
// File NoAdvection.cpp
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
// Description: Evaluation of the Navier Stokes advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/AdvectionTerms/NoAdvection.h>

using namespace std;

namespace Nektar
{

string NoAdvection::className
    = SolverUtils::GetAdvectionFactory().RegisterCreatorFunction(
            "NoAdvection",
            NoAdvection::create);

/**
 *
 */
NoAdvection::NoAdvection():
        Advection()

{
}


/**
 *
 */
NoAdvection::~NoAdvection()
{
}


/**
 *
 */
void NoAdvection::v_InitObject(
    LibUtilities::SessionReaderSharedPtr        pSession,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
{
}


/**
 *
 */
void NoAdvection::v_Advect(
    const int nConvectiveFields,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble> >        &advVel,
    const Array<OneD, Array<OneD, NekDouble> >        &inarray,
          Array<OneD, Array<OneD, NekDouble> >        &outarray,
    const NekDouble                                   &time,
    const Array<OneD, Array<OneD, NekDouble> >        &pFwd,
    const Array<OneD, Array<OneD, NekDouble> >        &pBwd)
{
    int nPointsTot = fields[0]->GetNpoints();
    for (int i = 0; i < inarray.size(); ++i)
    {
        Vmath::Zero(nPointsTot,outarray[i],1);
    }
}

} //end of namespace

