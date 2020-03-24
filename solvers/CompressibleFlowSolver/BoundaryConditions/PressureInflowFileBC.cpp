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

#include <boost/core/ignore_unused.hpp>

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
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    int nvariables = m_fields.size();
    // Loop over Boundary Regions for PressureInflowFileBC
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

void PressureInflowFileBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    boost::ignore_unused(time);

    int i, j;
    int nTracePts = m_fields[0]->GetTrace()->GetNpoints();
    int nVariables = physarray.size();
    int nDimensions = m_spacedim;

    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

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
    m_varConv->GetAbsoluteVelocity(Fwd, absVel);

    // Get speed of sound
    Array<OneD, NekDouble > pressure  (nTracePts);
    Array<OneD, NekDouble > soundSpeed(nTracePts);

    m_varConv->GetPressure(Fwd, pressure);
    m_varConv->GetSoundSpeed(Fwd, soundSpeed);

    // Get Mach
    Array<OneD, NekDouble > Mach(nTracePts, 0.0);
    Vmath::Vdiv(nTracePts, Vn, 1, soundSpeed, 1, Mach, 1);
    Vmath::Vabs(nTracePts, Mach, 1, Mach, 1);

    // Auxiliary variables
    int e, id1, id2, npts, pnt;
    NekDouble rhoeb;

    // Loop on the m_bcRegions
    for (e = 0; e < m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
         GetExpSize(); ++e)
    {
        npts = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
        GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
        GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        // Get internal energy
        Array<OneD, NekDouble> tmpPressure (npts, pressure+id2);
        Array<OneD, NekDouble> rho         (npts, Fwd[0]+id2);
        Array<OneD, NekDouble> e(npts);
        m_varConv->GetEFromRhoP(rho, tmpPressure, e);

        // Loop on points of m_bcRegion 'e'
        for (i = 0; i < npts; i++)
        {
            pnt = id2+i;

            // Subsonic flows
            if (Mach[pnt] < 0.99)
            {
                // Partial extrapolation for subsonic cases
                for (j = 0; j < nVariables-1; ++j)
                {
                    (m_fields[j]->GetBndCondExpansions()[m_bcRegion]->
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

                rhoeb = Fwd[0][pnt] * e[i] + Ek;

                (m_fields[nVariables-1]->GetBndCondExpansions()[m_bcRegion]->
                 UpdatePhys())[id1+i] = rhoeb;
            }
            // Supersonic flows
            else
            {
                for (j = 0; j < nVariables; ++j)
                {
                    // Extrapolation for supersonic cases
                    (m_fields[j]->GetBndCondExpansions()[m_bcRegion]->
                     UpdatePhys())[id1+i] = Fwd[j][pnt];
                }
            }
        }
    }
}

}
