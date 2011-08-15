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
