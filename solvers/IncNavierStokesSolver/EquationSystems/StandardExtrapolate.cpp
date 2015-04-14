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
     * 
     */
    void StandardExtrapolate::v_SubSteppingTimeIntegration(
        int intMethod,
        const LibUtilities::TimeIntegrationWrapperSharedPtr &IntegrationScheme)
    {
        switch(intMethod)
        {
            case LibUtilities::eIMEXOrder1:
            {
                m_intSteps = 1; 
            }
            break;
            case LibUtilities::eIMEXOrder2:
            {
                m_intSteps = 2;
            }
            break;
            case LibUtilities::eIMEXOrder3:
            {
                m_intSteps = 3;
            }
            break;
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
        const LibUtilities::TimeIntegrationSolutionSharedPtr &integrationSoln, 
        int nstep, 
        NekDouble time)
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

