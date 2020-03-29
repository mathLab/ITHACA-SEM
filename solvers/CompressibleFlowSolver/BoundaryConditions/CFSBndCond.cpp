///////////////////////////////////////////////////////////////////////////////
//
// File: CFSBndCond.cpp
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
// Description: Abstract base class for compressible solver boundary conditions.
//
///////////////////////////////////////////////////////////////////////////////

#include "CFSBndCond.h"

using namespace std;

namespace Nektar
{
CFSBndCondFactory& GetCFSBndCondFactory()
{
    static CFSBndCondFactory instance;
    return instance;
}

CFSBndCond::CFSBndCond(const LibUtilities::SessionReaderSharedPtr& pSession,
                const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
                const Array<OneD, Array<OneD, NekDouble> >&       pTraceNormals,
                const int pSpaceDim,
                const int bcRegion,
                const int cnt)
        : m_session(pSession),
        m_fields(pFields),
        m_traceNormals(pTraceNormals),
        m_spacedim(pSpaceDim),
        m_bcRegion(bcRegion),
        m_offset(cnt)
{
    m_velInf = Array<OneD, NekDouble> (m_spacedim, 0.0);
    m_session->LoadParameter("Gamma", m_gamma, 1.4);
    m_session->LoadParameter("rhoInf", m_rhoInf, 1.225);
    m_session->LoadParameter("pInf", m_pInf, 101325);
    m_session->LoadParameter("uInf", m_velInf[0], 0.1);
    if (m_spacedim >= 2)
    {
        m_session->LoadParameter("vInf", m_velInf[1], 0.0);
    }
    if (m_spacedim == 3)
    {
        m_session->LoadParameter("wInf", m_velInf[2], 0.0);
    }

    // Create auxiliary object to convert variables
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
                m_session, m_spacedim);
    
    m_diffusionAveWeight = 1.0;
}

/**
 * @param   bcRegion      id of the boundary region
 * @param   cnt           
 * @param   Fwd    
 * @param   physarray
 * @param   time
 */
void CFSBndCond::Apply(
            Array<OneD, Array<OneD, NekDouble> >               &Fwd,
            Array<OneD, Array<OneD, NekDouble> >               &physarray,
            const NekDouble                                    &time)
{
    v_Apply(Fwd, physarray, time);
}

/**
 * @ brief apply boundaryconditions for flow field derivatives
 * Currently only the ExtrapOrder0BC are used for all boundaries 
 * according to Cheng, Yang, Liu, JCP 327 (2016)484-502
 * @param   bcRegion      id of the boundary region
 * @param   cnt           
 * @param   Fwd    
 * @param   physarray
 * @param   time
 */
void CFSBndCond::ApplyDeriv(
    const Array<OneD, const Array<OneD, NekDouble> >               &Fwd,
    const Array<OneD, const Array<OneD, NekDouble> >               &physarray,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > > &DervFwd,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > > &dervarray,
    NekDouble                                                      time)
{
    v_ApplyDeriv(Fwd, physarray, DervFwd, dervarray, time);
}

void CFSBndCond::v_ApplyDeriv(
    const Array<OneD, const Array<OneD, NekDouble> >               &Fwd,
    const Array<OneD, const Array<OneD, NekDouble> >               &physarray,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > > &DervFwd,
    const Array<OneD, const Array<OneD, Array<OneD, NekDouble> > > &dervarray,
    NekDouble                                                      time)
{
    boost::ignore_unused(Fwd,dervarray,time);
    int i, j;
    int e, pnt;
    int id1, id2, nBCEdgePts;
    int nVariables = physarray.num_elements();
    int nDimensions = m_spacedim;

    const Array<OneD, const int> &traceBndMap
        = m_fields[0]->GetTraceBndMap();

    int eMax;

    eMax = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetExpSize();

    // Loop on m_bcRegions
    for (e = 0; e < eMax; ++e)
    {
        nBCEdgePts = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetExp(e)->GetTotPoints();
        id1 = m_fields[0]->GetBndCondExpansions()[m_bcRegion]->
            GetPhys_Offset(e) ;
        id2 = m_fields[0]->GetTrace()->GetPhys_Offset(traceBndMap[m_offset+e]);

        // Loop on points of m_bcRegion 'e'
        for (i = 0; i < nBCEdgePts; i++)
        {
            pnt = id2+i;

            // Setting up bcs for density and velocities
            for (int nd = 0; nd <nDimensions; ++nd)
            {
                for (j = 0; j <=nDimensions; ++j)
                {
                    (m_fields[j]->GetDerivBndCondExpansions()[m_bcRegion][nd]->
                    UpdatePhys())[id1+i] = DervFwd[nd][j][pnt];
                }

                // Setting up bcs for energy
                (m_fields[nVariables-1]->
                    GetDerivBndCondExpansions()[m_bcRegion][nd]->
                    UpdatePhys())[id1+i] = DervFwd[nd][nVariables-1][pnt];

            }
        }
    }
}

/**
 * @ brief Newly added bc should specify this virtual function
 * if the Bwd/value in m_bndCondExpansions is the target value like Direchlet bc weight should be 1.0.
 * if some average Fwd and Bwd/value in m_bndCondExpansions 
 * is the target value like WallViscousBC weight should be 0.5.
 */
void CFSBndCond::v_ApplyBwdWeight()
{
    NekDouble   weight  =   m_diffusionAveWeight;
    int nVariables = m_fields.num_elements();
    for(int i=0;i<nVariables;i++)
    {
        m_fields[i]->SetBndCondBwdWeight(m_bcRegion,weight);
    }
}

}
