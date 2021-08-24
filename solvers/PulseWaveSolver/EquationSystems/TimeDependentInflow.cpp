///////////////////////////////////////////////////////////////////////////////
//
// File CommMpi.cpp
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
// Description: TimeDependentInflow class
//
///////////////////////////////////////////////////////////////////////////////

#include <PulseWaveSolver/EquationSystems/TimeDependentInflow.h>

using namespace std;

namespace Nektar
{

std::string TimeDependentInflow::className =
    GetBoundaryFactory().RegisterCreatorFunction(
        "TimeDependent", TimeDependentInflow::create,
        "TimeDependent inflow boundary condition");

TimeDependentInflow::TimeDependentInflow(
    Array<OneD, MultiRegions::ExpListSharedPtr> pVessel,
    const LibUtilities::SessionReaderSharedPtr pSession,
    PulseWavePressureAreaSharedPtr pressureArea)
    : PulseWaveBoundary(pVessel, pSession, pressureArea)
{
}

TimeDependentInflow::~TimeDependentInflow()
{
}

void TimeDependentInflow::v_DoBoundary(
    const Array<OneD, const Array<OneD, NekDouble> > &inarray,
    Array<OneD, Array<OneD, NekDouble> > &A_0,
    Array<OneD, Array<OneD, NekDouble> > &beta,
    Array<OneD, Array<OneD, NekDouble> > &alpha,
    const NekDouble time, int omega, int offset, int n)
{
    Array<OneD, MultiRegions::ExpListSharedPtr> vessel(2);

    // Pointers to the domains
    vessel[0] = m_vessels[2 * omega];
    vessel[1] = m_vessels[2 * omega + 1];

    /* Evaluate the boundary conditions. Note that this does not assign new
       variables in the virtual region, so it only prescribes the inflow
       characteristic. */
    for (int i = 0; i < 2; ++i)
    {
        vessel[i]->EvaluateBoundaryConditions(time);
    }

}

} // namespace Nektar
