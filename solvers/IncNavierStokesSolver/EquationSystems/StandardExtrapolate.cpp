///////////////////////////////////////////////////////////////////////////////
//
// File: StandardExtrapolate.cpp
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
// Description: Abstract base class for StandardExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/StandardExtrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string StandardExtrapolate::className = GetExtrapolateFactory().RegisterCreatorFunction(
        "Standard",
        StandardExtrapolate::create,
        "Standard");

    StandardExtrapolate::StandardExtrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject)
        : Extrapolate(pSession,pFields,pPressure,pVel,advObject)
    {
    }

    StandardExtrapolate::~StandardExtrapolate()
    {
    }


    /** 
     * Function to extrapolate the new pressure boundary condition.
     * Based on the velocity field and on the advection term.
     * Acceleration term is also computed.
     * This routine is a general one for 2d and 3D application and it can be called
     * directly from velocity correction scheme. Specialisation on dimensionality is
     * redirected to the CalcNeumannPressureBCs method.
     */
    void StandardExtrapolate::v_EvaluatePressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        m_pressureCalls++;
        if(m_HBCnumber>0)
        {
            // Calculate non-linear and viscous BCs at current level
            // and put in m_pressureHBCs[0]
            CalcNeumannPressureBCs(fields,N,kinvis);
            
            // Extrapolate to n+1
            ExtrapolateArray(m_pressureHBCs);
            
            // Add (phi,Du/Dt) term to m_presureHBC
            AddDuDt();

            // Copy m_pressureHBCs to m_PbndExp
            CopyPressureHBCsToPbndExp();            
        }

        CalcOutflowBCs(fields, kinvis);
    }

    
    /** 
     * 
     */
    void StandardExtrapolate::v_SubSteppingTimeIntegration(
        const LibUtilities::TimeIntegrationSchemeSharedPtr & IntegrationScheme )
    {
        if ( IntegrationScheme->GetName() == "IMEX" ||
             IntegrationScheme->GetName() == "IMEXGear" )
        {
            m_intSteps = IntegrationScheme->GetOrder();
        }
        else
        {
            NEKERROR(ErrorUtil::efatal, "Integration method not suitable: "
                     "Options include IMEXGear or IMEXOrder{1,2,3,4}");
        }
    }

    /** 
     * 
     */
    void StandardExtrapolate::v_SubStepSetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &inarray, 
        NekDouble Aii_DT,
        NekDouble kinvis)
    {
    }


    /** 
     * 
     */
    void StandardExtrapolate::v_SubStepAdvance(
        const LibUtilities::TimeIntegrationScheme::TimeIntegrationSolutionSharedPtr & integrationSoln, 
              int                                                                     nstep, 
              NekDouble                                                               time )
    {
    }


    /** 
     * 
     */
    void StandardExtrapolate::v_SubStepSaveFields(
        int nstep)
    {
    }

    /** 
     * 
     */
    void StandardExtrapolate::v_MountHOPBCs(
        int HBCdata, 
        NekDouble kinvis, 
        Array<OneD, NekDouble> &Q, 
        Array<OneD, const NekDouble> &Advection)
    {
        Vmath::Svtvp(HBCdata,-kinvis,Q,1,Advection,1,Q,1);
    }
}
