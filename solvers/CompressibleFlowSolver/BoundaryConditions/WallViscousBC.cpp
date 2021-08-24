///////////////////////////////////////////////////////////////////////////////
//
// File: WallViscousBC.cpp
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
// Description: No-slip wall boundary condition
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include "WallViscousBC.h"

using namespace std;

namespace Nektar
{

std::string WallViscousBC::classNameViscous = GetCFSBndCondFactory().
    RegisterCreatorFunction("WallViscous",
                            WallViscousBC::create,
                            "No-slip (viscous) wall boundary condition.");

std::string WallViscousBC::classNameAdiabatic = GetCFSBndCondFactory().
    RegisterCreatorFunction("WallAdiabatic",
                            WallViscousBC::create,
                            "Adiabatic wall boundary condition.");

WallViscousBC::WallViscousBC(const LibUtilities::SessionReaderSharedPtr& pSession,
           const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
           const Array<OneD, Array<OneD, NekDouble> >& pTraceNormals,
           const int pSpaceDim,
           const int bcRegion,
           const int cnt)
    : CFSBndCond(pSession, pFields, pTraceNormals, pSpaceDim, bcRegion, cnt)
{
    m_diffusionAveWeight = 0.5;

    m_bndPhys = Array<OneD, Array<OneD, NekDouble> > (m_fields.size());
}

void WallViscousBC::v_Apply(
        Array<OneD, Array<OneD, NekDouble> >               &Fwd,
        Array<OneD, Array<OneD, NekDouble> >               &physarray,
        const NekDouble                                    &time)
{
    boost::ignore_unused(time);

    int i;
    int nVariables = physarray.size();

    // Find the fields whose WallViscous/Adiabatic-BC is time-dependent
    // Update variables on the boundaries of these fields
    // Get the updated variables on the WallViscous/Adiabatic boundary
    //
    // Maybe the EvaluateBoundaryConditions() should be put upstream to
    // CompressibleFlowSystem::NumCalRiemFluxJac(), So that the BCs will not
    // be repeatedly updated when there are more than one time-dependent BC.
    std::string varName;
    for (i = 0; i < nVariables; ++i)
    {
        if (m_fields[i]->GetBndConditions()[m_bcRegion]->IsTimeDependent())
        {
            varName = m_session->GetVariable(i);
            m_fields[i]->EvaluateBoundaryConditions(time, varName);

            m_bndPhys[i] = m_fields[i]->GetBndCondExpansions()[m_bcRegion]
                           ->UpdatePhys();
        }
    }


    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    // Take into account that for PDE based shock capturing, eps = 0 at the
    // wall. Adjust the physical values of the trace to take user defined
    // boundaries into account
    int e, id1, id2, nBCEdgePts, eMax;

    eMax = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize();

    for (e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetPhys_Offset(e);
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        // Boundary condition for epsilon term.
        if (nVariables == m_spacedim+3)
        {
            Vmath::Zero(nBCEdgePts, &Fwd[nVariables-1][id2], 1);
        }

        // V = - Vin
        for (i = 0; i < m_spacedim; i++)
        {
            Vmath::Neg(nBCEdgePts, &Fwd[i+1][id2], 1);
        }

        // Superimpose the perturbation
        for (i = 0; i < nVariables; ++i)
        {
            if (m_fields[i]->GetBndConditions()[m_bcRegion]->IsTimeDependent())
            {
                Vmath::Vadd(nBCEdgePts, &m_bndPhys[i][id1], 1,
                            &Fwd[i][id2], 1,  &Fwd[i][id2], 1);
            }
        }

        // Copy boundary adjusted values into the boundary expansion
        for (i = 0; i < nVariables; ++i)
        {
            Vmath::Vcopy(nBCEdgePts, &Fwd[i][id2], 1,
                         &(m_fields[i]->GetBndCondExpansions()[m_bcRegion]->
                           UpdatePhys())[id1], 1);
        }
    }
}

}
