///////////////////////////////////////////////////////////////////////////////
//
// File: AdvectionSystem.cpp
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
// Description: Base class for advection-based equation systems.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/AdvectionSystem.h>

namespace Nektar {
namespace SolverUtils {

/**
 *
 */
AdvectionSystem::AdvectionSystem(
    const LibUtilities::SessionReaderSharedPtr& pSession,
    const SpatialDomains::MeshGraphSharedPtr& pGraph)
    : UnsteadySystem(pSession, pGraph)
{
}


/**
 *
 */
AdvectionSystem::~AdvectionSystem()
{

}


/**
 *
 */
void AdvectionSystem::v_InitObject()
{
    UnsteadySystem::v_InitObject();
    m_session->LoadParameter("IO_CFLSteps", m_cflsteps, 0);
    m_session->LoadParameter("IO_CFLWriteFld", m_cflWriteFld, 0);
    m_session->LoadParameter("IO_CFLWriteFldWaitSteps", m_cflWriteFldWaitSteps, 0);
}

/**
 *
 */
bool AdvectionSystem::v_PostIntegrate(int step)
{
    bool result = UnsteadySystem::v_PostIntegrate(step);

    if((m_cflsteps && !((step+1)%m_cflsteps)) || m_cflWriteFld>0)
    {
        int elmtid;
        NekDouble cfl = GetCFLEstimate(elmtid);

        if(m_cflsteps && !((step+1)%m_cflsteps) && m_comm->GetRank() == 0)
        {
            if( m_HomogeneousType == eNotHomogeneous)
            {
                std::cout << "CFL: ";
            }
            else
            {
                std::cout << "CFL (zero plane): ";
            }
            std::cout << cfl << " (in elmt " << elmtid << ")" << std::endl;
        }
        
        // At each timestep, if cflWriteFld is set check if cfl is above treshold
        if(m_cflWriteFld>0 && cfl >= m_cflWriteFld && step >= m_cflWriteFldWaitSteps)
        {
            std::string outname =  m_sessionName +  "_CFLWriteFld";
            WriteFld(outname + ".fld");
            m_cflWriteFld = 0;            
        }
    }

    return result;
}

/**
 *
 */
Array<OneD, NekDouble>  AdvectionSystem::GetElmtCFLVals(void)
{
    int nelmt = m_fields[0]->GetExpSize();

    const Array<OneD, int> expOrder = GetNumExpModesPerExp();

    const NekDouble cLambda = 0.2; // Spencer book pag. 317

    Array<OneD, NekDouble> stdVelocity(nelmt, 0.0);
    stdVelocity = v_GetMaxStdVelocity();

    Array<OneD, NekDouble> cfl(nelmt, 0.0);
    NekDouble order;
    for(int el = 0; el < nelmt; ++el)
    {
        order = std::max(expOrder[el]-1, 1);
        cfl[el] =  m_timestep*(stdVelocity[el] * cLambda * order * order);
    }

    return cfl;
}

/**
 *
 */
NekDouble AdvectionSystem::GetCFLEstimate(int &elmtid)
{
    int n_element = m_fields[0]->GetExpSize();

    Array<OneD, NekDouble> cfl = GetElmtCFLVals();

    elmtid = Vmath::Imax(n_element,cfl,1);

    NekDouble CFL,CFL_loc;
    CFL = CFL_loc = cfl[elmtid];
    m_comm->AllReduce(CFL,LibUtilities::ReduceMax);

    // unshuffle elmt id if data is not stored in consecutive order.
    elmtid = m_fields[0]->GetExp(elmtid)->GetGeom()->GetGlobalID();
    if(CFL != CFL_loc)
    {
        elmtid = -1;
    }

    m_comm->AllReduce(elmtid,LibUtilities::ReduceMax);

    // express element id with respect to plane
    if(m_HomogeneousType == eHomogeneous1D)
    {
        elmtid = elmtid%m_fields[0]->GetPlane(0)->GetExpSize();
    }
    return CFL;
}

}
}
