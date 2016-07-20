///////////////////////////////////////////////////////////////////////////////
//
// File: PressureInflowFileBC.cpp
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
// Description: Pressure inflow boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include "PressureInflowFileBC.h"

using namespace std;

namespace Nektar
{

std::string PressureInflowFileBC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("PressureInflowFile",
                            PressureInflowFileBC::create,
                            "Pressure inflow (file) boundary condition.");

PressureInflowFileBC::PressureInflowFileBC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion)
{
    int nvariables = m_fields.num_elements();
    // Loop over Boundary Regions for PressureInflowFileBC
    m_fieldStorage = Array<OneD, Array<OneD, NekDouble> > (nvariables);

    int numBCPts = m_fields[0]->
        GetBndCondExpansions()[bcRegion]->GetNpoints();
    for (int i = 0; i < nvariables; ++i)
    {
        m_fieldStorage[i] = Array<OneD, NekDouble>(numBCPts, 0.0);
        Vmath::Vcopy(
            numBCPts,
            m_fields[i]->GetBndCondExpansions()[bcRegion]->GetPhys(), 1,
            m_fieldStorage[i], 1);
    }
}

void PressureInflowFileBC::v_Apply(
        int                                                 bcRegion,
        int                                                 cnt,
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    int i, j;
    int nTracePts = m_fields[0]->GetTrace()->GetNpoints();
    int nVariables = physarray.num_elements();
    int nDimensions = m_spacedim;

    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    NekDouble gammaMinusOne    = m_gamma - 1.0;
    NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

    Array<OneD, NekDouble> tmp1 (nTracePts, 0.0);
    Array<OneD, NekDouble> VnInf(nTracePts, 0.0);

    // Computing the normal velocity for characteristics coming
    // from outside the computational domain
    // Computing the normal velocity for characteristics coming
    // from outside the computational domain
    for( i =0; i < nDimensions; i++)
    {
        Vmath::Svtvp(nTracePts, m_velInf[i], 
                     m_traceNormals[i], 1,
                     VnInf, 1,
                     VnInf, 1);
    }

    // Computing the normal velocity for characteristics coming
    // from inside the computational domain
    Array<OneD, NekDouble > Vn (nTracePts, 0.0);
    Array<OneD, NekDouble > Vel(nTracePts, 0.0);
    for (i = 0; i < nDimensions; ++i)
    {
        Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, Vel, 1);
        Vmath::Vvtvp(nTracePts, m_traceNormals[i], 1, Vel, 1, Vn, 1, Vn, 1);
    }

    // Computing the absolute value of the velocity in order to compute the
    // Mach number to decide whether supersonic or subsonic
    Array<OneD, NekDouble > absVel(nTracePts, 0.0);
    for (i = 0; i < nDimensions; ++i)
    {
        Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp1, 1);
        Vmath::Vmul(nTracePts, tmp1, 1, tmp1, 1, tmp1, 1);
        Vmath::Vadd(nTracePts, tmp1, 1, absVel, 1, absVel, 1);
    }
    Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

    // Get speed of sound
    Array<OneD, NekDouble > soundSpeed (nTracePts, 0.0);
    Array<OneD, NekDouble > pressure   (nTracePts, 0.0);

    for (i = 0; i < nTracePts; i++)
    {
        if (m_spacedim == 1)
        {
            pressure[i] = (gammaMinusOne) * (Fwd[2][i] -
                            0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i]));
        }
        else if (m_spacedim == 2)
        {
            pressure[i] = (gammaMinusOne) * (Fwd[3][i] -
                            0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                   Fwd[2][i] * Fwd[2][i] / Fwd[0][i]));
        }
        else
        {
            pressure[i] = (gammaMinusOne) * (Fwd[4][i] -
                            0.5 * (Fwd[1][i] * Fwd[1][i] / Fwd[0][i] +
                                   Fwd[2][i] * Fwd[2][i] / Fwd[0][i] +
                                   Fwd[3][i] * Fwd[3][i] / Fwd[0][i]));
        }

        soundSpeed[i] = sqrt(m_gamma * pressure[i] / Fwd[0][i]);
    }

    // Get Mach
    Array<OneD, NekDouble > Mach(nTracePts, 0.0);
    Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
    Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);

    // Auxiliary variables
    int e, id1, id2, npts, pnt;
    NekDouble rhoeb;

    // Loop on the bcRegions
    for (e = 0; e < m_fields[0]->GetBndCondExpansions()[bcRegion]->
         GetExpSize(); ++e)
    {
        npts = m_fields[0]->GetBndCondExpansions()[bcRegion]->
        GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[bcRegion]->
        GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[cnt+e]);

        // Loop on points of bcRegion 'e'
        for (i = 0; i < npts; i++)
        {
            pnt = id2+i;

            // Subsonic flows
            if (Mach[pnt] < 0.99)
            {
                // Partial extrapolation for subsonic cases
                for (j = 0; j < nVariables-1; ++j)
                {
                    (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = m_fieldStorage[j][id1+i];
                }

                // Kinetic energy calculation
                NekDouble Ek = 0.0;
                for (j = 1; j < nVariables-1; ++j)
                {
                    Ek += 0.5 * (m_fieldStorage[j][id1+i] *
                                 m_fieldStorage[j][id1+i]) /
                                    m_fieldStorage[0][id1+i];
                }
                rhoeb = gammaMinusOneInv * pressure[pnt] + Ek;

                (m_fields[nVariables-1]->GetBndCondExpansions()[bcRegion]->
                 UpdatePhys())[id1+i] = rhoeb;
            }
            // Supersonic flows
            else
            {
                for (j = 0; j < nVariables; ++j)
                {
                    // Extrapolation for supersonic cases
                    (m_fields[j]->GetBndCondExpansions()[bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[j][pnt];
                }
            }
        }
    }
}

}
