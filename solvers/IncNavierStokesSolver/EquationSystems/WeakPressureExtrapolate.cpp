///////////////////////////////////////////////////////////////////////////////
//
// File: WeakPressureExtrapolate.cpp
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
// Description: Abstract base class for WeakPressureExtrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/WeakPressureExtrapolate.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
    /**
     * Registers the class with the Factory.
     */
    std::string WeakPressureExtrapolate::className = GetExtrapolateFactory().RegisterCreatorFunction(
        "WeakPressure",
        WeakPressureExtrapolate::create,
        "WeakPressure");

    WeakPressureExtrapolate::WeakPressureExtrapolate(
        const LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
        MultiRegions::ExpListSharedPtr pPressure,
        const Array<OneD, int> pVel,
        const SolverUtils::AdvectionSharedPtr advObject):
        StandardExtrapolate(pSession,pFields,pPressure,pVel,advObject)
    {
    }

    WeakPressureExtrapolate::~WeakPressureExtrapolate()
    {
    }


	
    /** 
     * Function to extrapolate the new pressure boundary condition.
     * Based on the velocity field and on the advection term.
     * Acceleration term is also computed.  This routine is a general
     * one for 2d and 3D application and it can be called directly
     * from velocity correction scheme. Specialisation on
     * dimensionality is redirected to the CalcNeumannPressureBCs
     * method.
     */
    void WeakPressureExtrapolate::v_EvaluatePressureBCs(
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        const Array<OneD, const Array<OneD, NekDouble> >  &N,
        NekDouble kinvis)
    {
        m_pressureCalls++;
        if(m_HBCnumber > 0)
        {
            // Calculate just viscous BCs at current level and put in
            // m_pressureHBCs[nlevels-1]
            CalcNeumannPressureBCs(fields,N,kinvis);
            
            // Extrapolate to m_pressureHBCs to n+1
            ExtrapolateArray(m_pressureHBCs);

            // \int_bnd q x n.u^{n} ds update current normal of field
            // add m_pressureHBCs to gamma_0/Dt * m_acceleration[0] 
            AddVelBC();
            
            // Copy m_pressureHBCs to m_PbndExp
            CopyPressureHBCsToPbndExp();            
        }
        // Evaluate High order outflow conditiosn if required. 
        CalcOutflowBCs(fields, kinvis);
    }

    void WeakPressureExtrapolate::AddVelBC(void)
    {
        int order = std::min(m_pressureCalls,m_intSteps);
        
        // Update velocity BF at n+1 (actually only needs doing if
        // velocity is time dependent on HBCs)
        IProductNormVelocityBCOnHBC(m_acceleration[0]);

        // Subtract acceleration term off m_pressureHBCs[nlevels-1]
        Vmath::Svtvp(m_numHBCDof,
                     -1.0*StifflyStable_Gamma0_Coeffs[order-1]/m_timestep,
                     m_acceleration[0],  1,
                     m_pressureHBCs[m_intSteps-1], 1,
                     m_pressureHBCs[m_intSteps-1], 1);
    }
	
    /** 
     *  vritual function which only puts in the curl operator into the bcs
     */
    void WeakPressureExtrapolate::v_MountHOPBCs(int HBCdata, 
                                                NekDouble kinvis, 
                                                Array<OneD, NekDouble> &Q, 
                                                Array<OneD, const NekDouble> &Advection)
    {
        Vmath::Smul(HBCdata,-kinvis,Q,1,Q,1);
    }
    

}

