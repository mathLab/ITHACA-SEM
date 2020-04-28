///////////////////////////////////////////////////////////////////////////////
//
// File: IsentropicVortexBC.cpp
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
// Description: Isentropic vortex boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "IsentropicVortexBC.h"

using namespace std;

namespace Nektar
{

std::string IsentropicVortexBC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("IsentropicVortex",
                            IsentropicVortexBC::create,
                            "Isentropic vortex boundary condition.");

IsentropicVortexBC::IsentropicVortexBC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
}

void IsentropicVortexBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    int nvariables = physarray.size();

    const Array<OneD, const int> &bndTraceMap = m_fields[0]->GetTraceBndMap();
    // loop over Boundary Regions
    int npoints, id1, id2, e_max;

    e_max = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize();

    for(int e = 0; e < e_max; ++e)
    {
        npoints = m_fields[0]->
            GetBndCondExpansions()[m_bcRegion]->GetExp(e)->GetTotPoints();
        id1  = m_fields[0]->
            GetBndCondExpansions()[m_bcRegion]->GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(bndTraceMap[m_offset+e]);

        Array<OneD,NekDouble> x(npoints, 0.0);
        Array<OneD,NekDouble> y(npoints, 0.0);
        Array<OneD,NekDouble> z(npoints, 0.0);

        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
        GetExp(e)->GetCoords(x, y, z);

        EvaluateIsentropicVortex(x, y, z, Fwd, time, id2);

        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Vcopy(npoints, &Fwd[i][id2], 1,
                         &(m_fields[i]->GetBndCondExpansions()[m_bcRegion]->
                         UpdatePhys())[id1], 1);
        }
    }
}

void IsentropicVortexBC::EvaluateIsentropicVortex(
    const Array<OneD, NekDouble>               &x,
    const Array<OneD, NekDouble>               &y,
    const Array<OneD, NekDouble>               &z,
    Array<OneD, Array<OneD, NekDouble> >       &u,
    NekDouble                                   time,
    const int                                   o)
{
    boost::ignore_unused(z);

    int nq = x.size();

    // Flow parameters
    const NekDouble x0    = 5.0;
    const NekDouble y0    = 0.0;
    const NekDouble beta  = 5.0;
    const NekDouble u0    = 1.0;
    const NekDouble v0    = 0.5;
    const NekDouble gamma = m_gamma;
    NekDouble r, xbar, ybar, tmp;
    NekDouble fac = 1.0/(16.0*gamma*M_PI*M_PI);

    // In 3D zero rhow field.
    if (m_spacedim == 3)
    {
        Vmath::Zero(nq, &u[3][o], 1);
    }

    // Fill storage
    for (int i = 0; i < nq; ++i)
    {
        xbar      = x[i] - u0*time - x0;
        ybar      = y[i] - v0*time - y0;
        r         = sqrt(xbar*xbar + ybar*ybar);
        tmp       = beta*exp(1-r*r);
        u[0][i+o] = pow(1.0 - (gamma-1.0)*tmp*tmp*fac, 1.0/(gamma-1.0));
        u[1][i+o] = u[0][i+o]*(u0 - tmp*ybar/(2*M_PI));
        u[2][i+o] = u[0][i+o]*(v0 + tmp*xbar/(2*M_PI));
        u[m_spacedim+1][i+o] = pow(u[0][i+o], gamma)/(gamma-1.0) +
        0.5*(u[1][i+o]*u[1][i+o] + u[2][i+o]*u[2][i+o]) / u[0][i+o];
    }
}

}
