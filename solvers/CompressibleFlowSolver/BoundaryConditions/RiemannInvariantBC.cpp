///////////////////////////////////////////////////////////////////////////////
//
// File: RiemannInvariantBC.cpp
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
// Description: Riemann invariant boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include "RiemannInvariantBC.h"

using namespace std;

namespace Nektar
{

std::string RiemannInvariantBC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("RiemannInvariant",
                            RiemannInvariantBC::create,
                            "Riemann invariant boundary condition.");

RiemannInvariantBC::RiemannInvariantBC(
           const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
}

void RiemannInvariantBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    int i, j;
    int nTracePts = m_fields[0]->GetTrace()->GetNpoints();
    int nDimensions = m_spacedim;

    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    NekDouble gammaInv         = 1.0 / m_gamma;
    NekDouble gammaMinusOne    = m_gamma - 1.0;
    NekDouble gammaMinusOneInv = 1.0 / gammaMinusOne;

    Array<OneD, NekDouble> VnInf(nTracePts, 0.0);

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
    Array<OneD, NekDouble> tmp (nTracePts, 0.0);
    Array<OneD, NekDouble > absVel(nTracePts, 0.0);
    for (i = 0; i < nDimensions; ++i)
    {
        Vmath::Vdiv(nTracePts, Fwd[i+1], 1, Fwd[0], 1, tmp, 1);
        Vmath::Vmul(nTracePts, tmp, 1, tmp, 1, tmp, 1);
        Vmath::Vadd(nTracePts, tmp, 1, absVel, 1, absVel, 1);
    }
    Vmath::Vsqrt(nTracePts, absVel, 1, absVel, 1);

    // Get speed of sound
    Array<OneD, NekDouble > soundSpeed(nTracePts);
    Array<OneD, NekDouble > pressure  (nTracePts);

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
    int eMax;
    int e, id1, id2, nBCEdgePts, pnt;
    NekDouble cPlus, rPlus, cMinus, rMinus, VDBC, VNBC;
    Array<OneD, NekDouble> velBC(nDimensions, 0.0);
    Array<OneD, NekDouble> rhoVelBC(nDimensions, 0.0);
    NekDouble rhoBC, EBC, cBC, sBC, pBC;

    eMax = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize();

    // Loop on m_bcRegions
    for (e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetExp(e)->GetTotPoints();

        id1 = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        // Loop on the points of the m_bcRegion
        for (i = 0; i < nBCEdgePts; i++)
        {
            pnt = id2+i;

            // Impose inflow Riemann invariant
            if (Vn[pnt] <= 0.0)
            {
                // Subsonic flows
                if (Mach[pnt] < 1.00)
                {
                    // + Characteristic from inside
                    cPlus = sqrt(m_gamma * pressure[pnt] / Fwd[0][pnt]);
                    rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                    // - Characteristic from boundary
                    cMinus = sqrt(m_gamma * m_pInf / m_rhoInf);
                    rMinus = VnInf[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                }
                else
                {
                    // + Characteristic from inside
                    cPlus = sqrt(m_gamma * m_pInf / m_rhoInf);
                    rPlus = VnInf[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                    // + Characteristic from inside
                    cMinus = sqrt(m_gamma * m_pInf / m_rhoInf);
                    rMinus = VnInf[pnt] - 2.0 * cPlus * gammaMinusOneInv;
                }

                // Riemann boundary variables
                VNBC = 0.5 * (rPlus + rMinus);
                cBC = 0.25 * gammaMinusOne * (rPlus - rMinus);
                VDBC = VNBC - VnInf[pnt];

                // Thermodynamic boundary variables
                sBC = m_pInf / (pow(m_rhoInf, m_gamma));
                rhoBC = pow((cBC * cBC) / (m_gamma * sBC), gammaMinusOneInv);
                pBC = rhoBC * cBC * cBC * gammaInv;

                // Kinetic energy initialiasation
                NekDouble EkBC = 0.0;

                // Boundary velocities
                for ( j = 0; j < nDimensions; ++j)
                {
                    velBC[j] = m_velInf[j] + VDBC * m_traceNormals[j][pnt];
                    rhoVelBC[j] = rhoBC * velBC[j];
                    EkBC += 0.5 * rhoBC * velBC[j]*velBC[j];
                }

                // Boundary energy
                EBC = pBC * gammaMinusOneInv + EkBC;

                // Imposing Riemann Invariant boundary conditions
                (m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
                 UpdatePhys())[id1+i] = rhoBC;
                for (j = 0; j < nDimensions; ++j)
                {
                    (m_fields[j+1]->GetBndCondExpansions()[m_bcRegion]->
                     UpdatePhys())[id1+i] = rhoVelBC[j];
                }
                (m_fields[nDimensions+1]->GetBndCondExpansions()[m_bcRegion]->
                 UpdatePhys())[id1+i] = EBC;

            }
            else // Impose outflow Riemann invariant
            {
                // Subsonic flows
                if (Mach[pnt] < 1.00)
                {
                    // + Characteristic from inside
                    cPlus = sqrt(m_gamma * pressure[pnt] / Fwd[0][pnt]);
                    rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                    // - Characteristic from boundary
                    cMinus = sqrt(m_gamma * m_pInf / m_rhoInf);
                    rMinus = VnInf[pnt] - 2.0 * cMinus * gammaMinusOneInv;
                }
                else
                {
                    // + Characteristic from inside
                    cPlus = sqrt(m_gamma * pressure[pnt] / Fwd[0][pnt]);
                    rPlus = Vn[pnt] + 2.0 * cPlus * gammaMinusOneInv;

                    // + Characteristic from inside
                    cMinus = sqrt(m_gamma * pressure[pnt] / Fwd[0][pnt]);
                    rMinus = Vn[pnt] - 2.0 * cPlus * gammaMinusOneInv;
                }

                // Riemann boundary variables
                VNBC = 0.5 * (rPlus + rMinus);
                cBC = 0.25 * gammaMinusOne * (rPlus - rMinus);
                VDBC = VNBC - Vn[pnt];

                // Thermodynamic boundary variables
                sBC = pressure[pnt] / (pow(Fwd[0][pnt], m_gamma));
                rhoBC = pow((cBC * cBC) / (m_gamma * sBC), gammaMinusOneInv);
                pBC = rhoBC * cBC * cBC * gammaInv;

                // Kinetic energy initialiasation
                NekDouble EkBC = 0.0;

                // Boundary velocities
                for ( j = 0; j < nDimensions; ++j)
                {
                    velBC[j] = Fwd[j+1][pnt] / Fwd[0][pnt] +
                                VDBC * m_traceNormals[j][pnt];
                    rhoVelBC[j] = rhoBC * velBC[j];
                    EkBC += 0.5 * rhoBC * velBC[j]*velBC[j];
                }

                // Boundary energy
                EBC = pBC * gammaMinusOneInv + EkBC;

                // Imposing Riemann Invariant boundary conditions
                (m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
                 UpdatePhys())[id1+i] = rhoBC;
                for (j = 0; j < nDimensions; ++j)
                {
                    (m_fields[j+1]->GetBndCondExpansions()[m_bcRegion]->
                     UpdatePhys())[id1+i] = rhoVelBC[j];
                }
                (m_fields[nDimensions+1]->GetBndCondExpansions()[m_bcRegion]->
                 UpdatePhys())[id1+i] = EBC;
            }
        }
    }
}

}
