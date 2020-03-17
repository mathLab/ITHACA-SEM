///////////////////////////////////////////////////////////////////////////////
//
// File: SubSteppingExtrapolateWeakPressure.cpp
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
// Description: Abstract base class for SubSteppingExtrapolate with Weak Pressure VCS
//              mainly redefines SubStepSetPressureBCs
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/SubSteppingExtrapolateWeakPressure.h>
#include <LibUtilities/Communication/Comm.h>

using namespace std;

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string SubSteppingExtrapolateWeakPressure::className = GetExtrapolateFactory().RegisterCreatorFunction(
        "SubSteppingWeakPressure",
        SubSteppingExtrapolateWeakPressure::create,
        "SubSteppingWeakPressure");

    SubSteppingExtrapolateWeakPressure::SubSteppingExtrapolateWeakPressure(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr  pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject)
        : SubSteppingExtrapolate(pSession,pFields,pPressure,pVel,advObject)
    {
    }

    SubSteppingExtrapolateWeakPressure::~SubSteppingExtrapolateWeakPressure()
    {
    }

    void SubSteppingExtrapolateWeakPressure::v_SubStepSetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        const NekDouble Aii_Dt,
        NekDouble kinvis)
    {
        Array<OneD, Array<OneD, NekDouble> > nullvelfields;

        m_pressureCalls++;

        // Calculate non-linear and viscous BCs at current level and
        // put in m_pressureHBCs[0]
        CalcNeumannPressureBCs(inarray,nullvelfields,kinvis);

        // Extrapolate to m_pressureHBCs to n+1
        ExtrapolateArray(m_pressureHBCs);

        // Add (phi,gamma0 u^{n+1}/Dt) term to m_presureHBC
        AddVelBC();

        // Copy m_pressureHBCs to m_PbndExp
        CopyPressureHBCsToPbndExp();

        // Evaluate High order outflow conditiosn if required.
        CalcOutflowBCs(inarray, kinvis);
    }


    // In weak pressure formulation  we also require \int q u.n ds on outflow boundary
    void SubSteppingExtrapolateWeakPressure::v_AddNormVelOnOBC(const int cnt, const int nreg,
                                                               Array<OneD, Array<OneD, NekDouble> > &u)
    {
        if(!m_houtflow.get()) // no outflow on partition so just return
        {
           return;
        }

        int nbcoeffs = m_PBndExp[nreg]->GetNcoeffs();
        //int nqb      = m_PBndExp[nreg]->GetTotPoints();

        Array<OneD, NekDouble> IProdVnTmp(nbcoeffs);

        Array<OneD, Array<OneD, NekDouble> > ubnd(m_curl_dim);


        for(int i = 0; i < m_curl_dim; ++i)
        {
            EvaluateBDFArray(m_houtflow->m_outflowVelBnd[cnt][i]);

            ubnd[i] = m_houtflow->m_outflowVelBnd[cnt][i][m_intSteps-1];

            // point input u to the first part of the array for later uee.
            u[i] = m_houtflow->m_outflowVelBnd[cnt][i][0];
        }

        m_PBndExp[nreg]->NormVectorIProductWRTBase(ubnd,IProdVnTmp);

        Vmath::Svtvp(nbcoeffs,-1.0/m_timestep,IProdVnTmp,1,m_PBndExp[nreg]->UpdateCoeffs(),1,
                     m_PBndExp[nreg]->UpdateCoeffs(),1);
    }
}
