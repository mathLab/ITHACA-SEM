///////////////////////////////////////////////////////////////////////////////
//
// File DriverArnoldi.cpp
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
// Description: Base Driver class for the stability solver
//
///////////////////////////////////////////////////////////////////////////////

#include <Auxiliary/DriverArnoldi.h>

namespace Nektar
{
    /**
     *
     */
    DriverArnoldi::DriverArnoldi(LibUtilities::SessionReaderSharedPtr pSession)
            : Driver(pSession)
    {
        m_session->LoadParameter("IO_InfoSteps", m_infosteps, 0);
    };


    /**
     *
     */
    DriverArnoldi::~DriverArnoldi()
    {
    };


    /**
     *
     */
    void DriverArnoldi::v_InitObject()
    {
    	Driver::v_InitObject();
    	
        m_session->MatchSolverInfo("SolverType","VelocityCorrectionScheme",m_TimeSteppingAlgorithm, false);

        if(m_TimeSteppingAlgorithm)
        {
            m_period  = m_session->GetParameter("TimeStep")* m_session->GetParameter("NumSteps");
            m_nfields = m_equ[0]->UpdateFields().num_elements() - 1;
        }
        else
        {
            m_period  = 1.0;
            ASSERTL0(m_session->DefinesFunction("BodyForce"),"A BodyForce section needs to be defined for this solver type");
            m_nfields = m_equ[0]->UpdateFields().num_elements();
        }

        m_session->LoadParameter("kdim",  m_kdim,  16);
        m_session->LoadParameter("nvec",  m_nvec,  2);
        m_session->LoadParameter("nits",  m_nits,  500);
        m_session->LoadParameter("evtol", m_evtol, 1e-06);
    }


    /**
     * Copy Arnoldi array to field variables which depend from
     * either the m_fields or m_forces
     */
    void DriverArnoldi::CopyArnoldiArrayToField(Array<OneD, NekDouble> &array)
    {
        Array<OneD, MultiRegions::ExpListSharedPtr> fields;
        if(m_TimeSteppingAlgorithm)
        {
            fields = m_equ[0]->UpdateFields();
        }
        else
        {
            fields = m_equ[0]->UpdateForces();
        }

        int nq = fields[0]->GetNpoints();

        for (int k = 0; k < m_nfields; ++k)
        {
            Vmath::Vcopy(nq, &array[k*nq], 1, &fields[k]->UpdatePhys()[0], 1);
            fields[k]->SetPhysState(true);
        }
    };


    /**
     * Copy field variables which depend from either the m_fields
     * or m_forces array the Arnoldi array
     */
    void DriverArnoldi::CopyFieldToArnoldiArray(Array<OneD, NekDouble> &array)
    {
        Array<OneD, MultiRegions::ExpListSharedPtr> fields;
        fields = m_equ[m_nequ-1]->UpdateFields();
        int nq = fields[0]->GetNpoints();

        for (int k = 0; k < m_nfields; ++k)
        {
            if(!m_TimeSteppingAlgorithm)
            {
                fields[k]->BwdTrans_IterPerExp(fields[k]->GetCoeffs(),fields[k]->UpdatePhys());
            }

            Vmath::Vcopy(nq,  &fields[k]->GetPhys()[0], 1, &array[k*nq], 1);
            fields[k]->SetPhysState(true);
        }
    };
}
