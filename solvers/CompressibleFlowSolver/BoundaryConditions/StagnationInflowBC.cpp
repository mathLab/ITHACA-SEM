///////////////////////////////////////////////////////////////////////////////
//
// File: StagnationInflowBC.cpp
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
// Description: Stagnation conditions inflow boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "StagnationInflowBC.h"

using namespace std;

namespace Nektar
{

std::string StagnationInflowBC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("StagnationInflow",
                            StagnationInflowBC::create,
                            "Stagnation conditions inflow boundary condition.");

StagnationInflowBC::StagnationInflowBC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    int nvariables = m_fields.size();
    int expdim     = pFields[0]->GetGraph()->GetMeshDimension();
    int spacedim   = pFields[0]->GetGraph()->GetSpaceDimension();
    m_swirl        = ((spacedim == 3) && (expdim == 2));
    // Loop over Boundary Regions for StagnationInflowBC
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

void StagnationInflowBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    boost::ignore_unused(time);

    int i, j;
    int nTracePts  = m_fields[0]->GetTrace()->GetNpoints();
    int numBCPts   = m_fields[0]->
        GetBndCondExpansions()[m_bcRegion]->GetNpoints();
    int nVariables = physarray.size();

    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    NekDouble gammaInv = 1.0 / m_gamma;
    NekDouble gammaMinusOne    = m_gamma - 1.0;
    NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

    // Get stagnation pressure (with zero swirl)
    Array<OneD, NekDouble > pStag      (numBCPts);
    if (m_swirl)
    {
        Array<OneD, NekDouble > tmp       (numBCPts);
        Vmath::Vcopy(numBCPts, m_fieldStorage[3], 1, tmp, 1);
        Vmath::Zero(numBCPts, m_fieldStorage[3], 1);
        m_varConv->GetPressure(m_fieldStorage, pStag);
        Vmath::Vcopy(numBCPts, tmp, 1, m_fieldStorage[3], 1);
    }
    else
    {
        m_varConv->GetPressure(m_fieldStorage, pStag);
    }

    // Get Mach from Fwd
    Array<OneD, NekDouble > soundSpeed(nTracePts);
    Array<OneD, NekDouble > Mach      (nTracePts);
    m_varConv->GetSoundSpeed(Fwd, soundSpeed);
    m_varConv->GetMach(Fwd, soundSpeed, Mach);

    // Auxiliary variables
    int e, id1, id2, npts, pnt;
    NekDouble rho, p;

    // Loop on the m_bcRegion
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

            // Pressure from stagnation pressure and Mach
            p = pStag[id1+i] /
                pow(1.0 + (gammaMinusOne)/2.0 * Mach[pnt]*Mach[pnt],
                    m_gamma/(gammaMinusOne));

            // rho from isentropic relation: rho = rhoStag *(p/pStag)^1/Gamma
            rho = m_fieldStorage[0][id1+i] *
                    pow(p/pStag[id1+i],gammaInv);
            (m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
                UpdatePhys())[id1+i] = rho;

            // Extrapolation for velocity and Kinetic energy calculation
            int lim = m_swirl ? nVariables - 2 : nVariables - 1;
            NekDouble Ek = 0.0;
            for (j = 1; j < lim; ++j)
            {
                (m_fields[j]->GetBndCondExpansions()[m_bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[j][pnt];

                Ek += 0.5 * (Fwd[j][pnt] * Fwd[j][pnt]) / Fwd[0][pnt];
            }

            if (m_swirl)
            {
                // Prescribed swirl
                (m_fields[3]->GetBndCondExpansions()[m_bcRegion]->
                    UpdatePhys())[id1+i] = m_fieldStorage[3][id1+i];
                Ek += 0.5 * (m_fieldStorage[3][id1+i] *
                             m_fieldStorage[3][id1+i]) / Fwd[0][pnt];
            }

            // Energy
            (m_fields[nVariables-1]->GetBndCondExpansions()[m_bcRegion]->
                UpdatePhys())[id1+i] = p * gammaMinusOneInv + Ek;
        }
    }
}

}
