///////////////////////////////////////////////////////////////////////////////
//
// File: StagnationInflowSwirlBC.cpp
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
// Description: Stagnation conditions inflow with swirl boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include "StagnationInflowSwirlBC.h"

using namespace std;

namespace Nektar
{

std::string StagnationInflowSwirlBC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("StagnationInflowSwirl",
                            StagnationInflowSwirlBC::create,
                            "Stagnation conditions inflow with Swirl boundary condition.");

StagnationInflowSwirlBC::StagnationInflowSwirlBC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    int nvariables = m_fields.num_elements();
    // Loop over Boundary Regions for StagnationInflowSwirlBC
    m_fieldStorage = Array<OneD, Array<OneD, NekDouble> > (nvariables);

    int numBCPts = m_fields[0]->
        GetBndCondExpansions()[m_bcRegion]->GetNpoints();
    for (int i = 0; i < nvariables; ++i)
    {
        m_fieldStorage[i] = Array<OneD, NekDouble>(numBCPts, 0.0);
        Vmath::Vcopy(
            numBCPts,
            m_fields[i]->GetBndCondExpansions()[m_bcRegion]->GetPhys(), 1,
            m_fieldStorage[i], 1);
    }
}

void StagnationInflowSwirlBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    int i, j;
    int nTracePts = m_fields[0]->GetTrace()->GetNpoints();
    int numBCPts = m_fields[0]->
        GetBndCondExpansions()[m_bcRegion]->GetNpoints();
    int nVariables = physarray.num_elements();

    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    NekDouble gammaInv = 1.0 / m_gamma;

    // Get pressure
    Array<OneD, NekDouble > pressure  (nTracePts);
    m_varConv->GetPressure(Fwd, pressure);

    // Get reference pressure (with zero swirl)
    Array<OneD, NekDouble > pRef      (numBCPts);
    Array<OneD, NekDouble > tmp       (numBCPts);
    Vmath::Vcopy(numBCPts, m_fieldStorage[3], 1, tmp, 1);
    Vmath::Zero(numBCPts, m_fieldStorage[3], 1);
    m_varConv->GetPressure(m_fieldStorage, pRef);
    Vmath::Vcopy(numBCPts, tmp, 1, m_fieldStorage[3], 1);

    // Auxiliary variables
    int e, id1, id2, npts, pnt;
    NekDouble rho;

    // Loop on the m_bcRegions
    for (e = 0; e < m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
         GetExpSize(); ++e)
    {
        npts = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
        GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        // Loop on points of m_bcRegion 'e'
        for (i = 0; i < npts; i++)
        {
            pnt = id2 + i;
            // Density from isentropic relation: rho = rhoRef *(p/pRef)^1/Gamma
            rho = m_fieldStorage[0][id1+i] * 
                    pow(pressure[pnt]/pRef[id1+i],gammaInv);
            (m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
                UpdatePhys())[id1+i] = rho;

            // Extrapolation for in plane velocity
            for (j = 1; j < nVariables-2; ++j)
            {
                (m_fields[j]->GetBndCondExpansions()[m_bcRegion]->
                 UpdatePhys())[id1+i] = Fwd[j][pnt];
            }

            // Prescribed swirl
            (m_fields[nVariables-2]->GetBndCondExpansions()[m_bcRegion]->
                UpdatePhys())[id1+i] = m_fieldStorage[nVariables-2][id1+i];

            // Prescribed energy
            (m_fields[nVariables-1]->GetBndCondExpansions()[m_bcRegion]->
                UpdatePhys())[id1+i] = m_fieldStorage[nVariables-1][id1+i];
        }
    }
}

}
