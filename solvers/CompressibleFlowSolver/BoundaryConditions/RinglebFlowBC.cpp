///////////////////////////////////////////////////////////////////////////////
//
// File: RinglebFlowBC.cpp
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
// Description: Ringleb flow boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include "RinglebFlowBC.h"
#include <MultiRegions/AssemblyMap/AssemblyMapDG.h>

using namespace std;

namespace Nektar
{

std::string RinglebFlowBC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("RinglebFlow",
                            RinglebFlowBC::create,
                            "Ringleb flow boundary condition.");

RinglebFlowBC::RinglebFlowBC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    m_expdim = pFields[0]->GetGraph()->GetMeshDimension();

    m_homo1D = false;
    if (m_session->DefinesSolverInfo("HOMOGENEOUS"))
    {
        std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
        if ((HomoStr == "HOMOGENEOUS1D") || (HomoStr == "Homogeneous1D")
                        || (HomoStr == "1D") || (HomoStr == "Homo1D"))
        {
            m_homo1D = true;
        }
    }
}

void RinglebFlowBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    int nvariables      = physarray.size();

    // For 3DHomogenoeus1D
    int n_planes = 1;
    if (m_expdim == 2 &&  m_homo1D)
    {
        int nPointsTot = m_fields[0]->GetTotPoints();
        int nPointsTot_plane = m_fields[0]->GetPlane(0)->GetTotPoints();
        n_planes = nPointsTot/nPointsTot_plane;
    }

    int id2, id2_plane, e_max;

    e_max = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize();

    for(int e = 0; e < e_max; ++e)
    {
        int npoints = m_fields[0]->
            GetBndCondExpansions()[m_bcRegion]->GetExp(e)->GetTotPoints();
        int id1  = m_fields[0]->
            GetBndCondExpansions()[m_bcRegion]->GetPhys_Offset(e);

        // For 3DHomogenoeus1D
        if (m_expdim == 2 &&  m_homo1D)
        {
            int m_offset_plane = m_offset/n_planes;
            int e_plane;
            int e_max_plane = e_max/n_planes;
            int nTracePts_plane = m_fields[0]->GetTrace()->GetNpoints();

            int planeID = floor((e + 0.5 )/ e_max_plane );
            e_plane = e - e_max_plane*planeID;

            id2_plane  = m_fields[0]->GetTrace()->GetPhys_Offset(
                m_fields[0]->GetTraceMap()->
                GetBndCondIDToGlobalTraceID(m_offset_plane + e_plane));
            id2 = id2_plane + planeID*nTracePts_plane;
        }
        else // For general case
        {
            id2 = m_fields[0]->
                GetTrace()->GetPhys_Offset(m_fields[0]->GetTraceMap()->
                GetBndCondIDToGlobalTraceID(m_offset+e));
        }

        Array<OneD,NekDouble> x0(npoints, 0.0);
        Array<OneD,NekDouble> x1(npoints, 0.0);
        Array<OneD,NekDouble> x2(npoints, 0.0);

        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
        GetExp(e)->GetCoords(x0, x1, x2);

        // Flow parameters
        NekDouble c, k, phi, r, J, VV, pp, sint, P, ss;
        NekDouble J11, J12, J21, J22, det;
        NekDouble Fx, Fy;
        NekDouble xi, yi;
        NekDouble dV;
        NekDouble dtheta;
        NekDouble par1;
        NekDouble theta     = M_PI / 4.0;
        NekDouble kExt      = 0.7;
        NekDouble V         = kExt * sin(theta);
        NekDouble toll      = 1.0e-8;
        NekDouble errV      = 1.0;
        NekDouble errTheta  = 1.0;
        NekDouble gamma     = m_gamma;
        NekDouble gamma_1_2 = (gamma - 1.0) / 2.0;

        // Loop on all the points of that edge
        for (int j = 0; j < npoints; j++)
        {

            while ((abs(errV) > toll) || (abs(errTheta) > toll))
            {
                VV   = V * V;
                sint = sin(theta);
                c    = sqrt(1.0 - gamma_1_2 * VV);
                k    = V / sint;
                phi  = 1.0 / k;
                pp   = phi * phi;
                J    = 1.0 / c + 1.0 / (3.0 * c * c * c) +
                1.0 / (5.0 * c * c * c * c * c) -
                0.5 * log((1.0 + c) / (1.0 - c));

                r    = pow(c, 1.0 / gamma_1_2);
                xi   = 1.0 / (2.0 * r) * (1.0 / VV - 2.0 * pp) + J / 2.0;
                yi   = phi / (r * V) * sqrt(1.0 - VV * pp);
                par1 = 25.0 - 5.0 * VV;
                ss   = sint * sint;

                Fx   = xi - x0[j];
                Fy   = yi - x1[j];

                J11 = 39062.5 / pow(par1, 3.5) *
                (1.0 / VV - 2.0 / VV * ss) * V + 1562.5 /
                pow(par1, 2.5) * (-2.0 / (VV * V) + 4.0 /
                (VV * V) * ss) + 12.5 / pow(par1, 1.5) * V +
                312.5 / pow(par1, 2.5) * V + 7812.5 /
                pow(par1, 3.5) * V - 0.25 *
                (-1.0 / pow(par1, 0.5) * V / (1.0 - 0.2 *
                pow(par1, 0.5)) - (1.0 + 0.2 * pow(par1, 0.5)) /
                pow((1.0 - 0.2 * pow(par1, 0.5)), 2.0) /
                pow(par1, 0.5) * V) / (1.0 + 0.2 * pow(par1, 0.5)) *
                (1.0 - 0.2 * pow(par1, 0.5));

                J12 = -6250.0 / pow(par1, 2.5) / VV * sint * cos(theta);
                J21 = -6250.0 / (VV * V) * sint / pow(par1, 2.5) *
                pow((1.0 - ss), 0.5) + 78125.0 / V * sint /
                pow(par1, 3.5) * pow((1.0 - ss), 0.5);

                // the matrix is singular when theta = pi/2
                if (abs(x1[j]) < toll && abs(cos(theta)) < toll)
                {
                    J22 = -39062.5 / pow(par1, 3.5) / V + 3125 /
                    pow(par1, 2.5) / (VV * V) + 12.5 /
                    pow(par1, 1.5) * V + 312.5 / pow(par1, 2.5) *
                    V + 7812.5 / pow(par1, 3.5) * V - 0.25 *
                    (-1.0 / pow(par1, 0.5) * V / (1.0 - 0.2 *
                    pow(par1, 0.5)) - (1.0 + 0.2 * pow(par1, 0.5)) /
                    pow((1.0 - 0.2 * pow(par1, 0.5)), 2.0) /
                    pow(par1, 0.5) * V) / (1.0 + 0.2 *
                    pow(par1, 0.5)) * (1.0 - 0.2 * pow(par1, 0.5));

                    // dV = -dV/dx * Fx
                    dV      = -1.0 / J22 * Fx;
                    dtheta  = 0.0;
                    theta   = M_PI / 2.0;
                }
                else
                {
                    J22 = 3125.0 / VV * cos(theta) / pow(par1, 2.5) *
                    pow((1.0 - ss), 0.5) - 3125.0 / VV * ss /
                    pow(par1, 2.5) / pow((1.0 - ss), 0.5) *
                    cos(theta);

                    det = -1.0 / (J11 * J22 - J12 * J21);

                    // [dV dtheta]' = -[invJ]*[Fx Fy]'
                    dV     = det * ( J22 * Fx - J12 * Fy);
                    dtheta = det * (-J21 * Fx + J11 * Fy);
                }

                V     = V + dV;
                theta = theta + dtheta;

                errV     = abs(dV);
                errTheta = abs(dtheta);
            }

            c                      = sqrt(1.0 - gamma_1_2 * VV);
            int kk                 = id2 + j;
            NekDouble timeramp     = 200.0;;
            if (time<timeramp &&
                !(m_session->DefinesFunction("InitialConditions") &&
                m_session->GetFunctionType("InitialConditions", 0) ==
                LibUtilities::eFunctionTypeFile))
            {
                Fwd[0][kk] = pow(c, 1.0 / gamma_1_2) *
                exp(-1.0 + time /timeramp);

                Fwd[1][kk] = Fwd[0][kk] * V * cos(theta) *
                exp(-1 + time / timeramp);

                Fwd[2][kk] = Fwd[0][kk] * V * sin(theta) *
                exp(-1 + time / timeramp);
            }
            else
            {
                Fwd[0][kk]  = pow(c, 1.0 / gamma_1_2);
                Fwd[1][kk] = Fwd[0][kk] * V * cos(theta);
                Fwd[2][kk] = Fwd[0][kk] * V * sin(theta);
            }

            P  = (c * c) * Fwd[0][kk] / gamma;
            Fwd[3][kk]  = P / (gamma - 1.0) + 0.5 *
            (Fwd[1][kk] * Fwd[1][kk] / Fwd[0][kk] +
            Fwd[2][kk] * Fwd[2][kk] / Fwd[0][kk]);

            errV     = 1.0;
            errTheta = 1.0;
            theta    = M_PI / 4.0;
            V        = kExt * sin(theta);
        }

        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npoints, &Fwd[i][id2], 1,
                         &(m_fields[i]->GetBndCondExpansions()[m_bcRegion]->
                         UpdatePhys())[id1],1);
        }
    }
}

}
