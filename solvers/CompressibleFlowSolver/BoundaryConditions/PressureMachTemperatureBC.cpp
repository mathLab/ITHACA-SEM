///////////////////////////////////////////////////////////////////////////////
//
// File: PressureMachTemperatureBC.cpp
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
// Description: Boundary condition specified in terms of pressure, Mach number
//              and temperature
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "PressureMachTemperatureBC.h"

using namespace std;

namespace Nektar
{

std::string PressureMachTemperatureBC::className = GetCFSBndCondFactory().
    RegisterCreatorFunction("PressureMachTemperature",
                            PressureMachTemperatureBC::create,
                            "BC prescribed in terms of p, Ma and T.");

PressureMachTemperatureBC::PressureMachTemperatureBC(
           const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    int nvariables = m_fields.size();
    int numBCPts = m_fields[0]->
        GetBndCondExpansions()[m_bcRegion]->GetNpoints();

    // Array for storing conserved variables on the boundary
    m_bcStorage = Array<OneD, Array<OneD, NekDouble> > (nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        m_bcStorage[i] = Array<OneD, NekDouble> (numBCPts, 0.0);
    }

    // We assume that the pressure is given in entry [0] of
    //   the BC ("rho" position) and the temperature in entry m_spacedim+1
    //   ("E" position)
    const Array<OneD, const NekDouble> pressure =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetPhys();
    const Array<OneD, const NekDouble> temperature =
        m_fields[m_spacedim+1]->GetBndCondExpansions()[m_bcRegion]->GetPhys();

    // Calculate density
    m_varConv->GetRhoFromPT(pressure, temperature, m_bcStorage[0]);
    // Calculate the internal energy times density
    m_varConv->GetEFromRhoP(m_bcStorage[0], pressure,
                            m_bcStorage[m_spacedim+1]);
    Vmath::Vmul(numBCPts, m_bcStorage[m_spacedim+1], 1,
                          m_bcStorage[0], 1,
                          m_bcStorage[m_spacedim+1], 1);
    // We can now obtain the sound speed at this (rho,e) condition
    Array<OneD, NekDouble> soundSpeed (numBCPts);
    m_varConv->GetSoundSpeed(m_bcStorage, soundSpeed);

    // Now update momentum and add kinetic energy to E
    Array<OneD, NekDouble> tmp (numBCPts);
    for (int i = 0; i < m_spacedim; ++i)
    {
        // tmp = velocity in i direction
        Vmath::Vmul(numBCPts,
                m_fields[i+1]->GetBndCondExpansions()[m_bcRegion]->GetPhys(), 1,
                soundSpeed, 1,
                tmp, 1);
        // rho*u
        Vmath::Vmul(numBCPts,
                m_bcStorage[0], 1,
                tmp, 1,
                m_bcStorage[i+1], 1);
        // tmp = 0.5 * rho *(rhou) in vel
        Vmath::Vmul(numBCPts,
                m_bcStorage[i+1], 1,
                tmp, 1,
                tmp, 1);
        Vmath::Smul(numBCPts,
                0.5,
                tmp, 1,
                tmp, 1);
        // Add to E
        Vmath::Vadd(numBCPts,
                m_bcStorage[m_spacedim+1], 1,
                tmp, 1,
                m_bcStorage[m_spacedim+1], 1);
    }

    // Copy to boundary condition
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(
            numBCPts,
            m_bcStorage[i], 1,
            m_fields[i]->GetBndCondExpansions()[m_bcRegion]->UpdatePhys(), 1);
    }
}

void PressureMachTemperatureBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    boost::ignore_unused(Fwd, physarray, time);

    int nvariables = m_fields.size();
    int numBCPts = m_fields[0]->
        GetBndCondExpansions()[m_bcRegion]->GetNpoints();
    // Copy conserved variables to boundary condition
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vcopy(
            numBCPts,
            m_bcStorage[i], 1,
            m_fields[i]->GetBndCondExpansions()[m_bcRegion]->UpdatePhys(), 1);
    }
}

}
