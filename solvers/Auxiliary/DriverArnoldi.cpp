/*
 * DriverArnoldi.cpp
 *
 *  Created on: 14 Aug 2011
 *      Author: cc
 */

#include <Auxiliary/DriverArnoldi.h>

namespace Nektar
{
    /**
     *
     */
    DriverArnoldi::DriverArnoldi(LibUtilities::SessionReaderSharedPtr pSession)
            : Driver(pSession)
    {
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
            //m_forces  = m_equ[0]->UpdateForces();
            m_nfields = m_equ[0]->UpdateFields().num_elements();
        }

        m_session->LoadParameter("kdim",  m_kdim,  8);
        m_session->LoadParameter("nvec",  m_nvec,  1);
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
