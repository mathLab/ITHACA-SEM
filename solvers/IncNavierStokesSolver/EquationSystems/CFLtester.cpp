///////////////////////////////////////////////////////////////////////////////
//
// File CLFtester.cpp
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
// Description: CFL tester solve routines
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

namespace Nektar
{
	
    NekDouble IncNavierStokes::v_GetTimeStep(const Array<OneD,int> ExpOrder, const Array<OneD,NekDouble> CFL, NekDouble timeCFL)
    { 
        
        int nvariables = m_fields.num_elements();   // Number of variables in the mesh
        int nTotQuadPoints  = GetTotPoints();
        int n_element  = m_fields[0]->GetExpSize(); 
        Array<OneD, NekDouble> tstep(n_element,0.0);
        const NekDouble minLengthStdTri  = 1.414213;
        const NekDouble minLengthStdQuad = 2.0;
        const NekDouble cLambda = 0.2; // Spencer book pag. 317
        Array<OneD, NekDouble> stdVelocity(n_element,0.0);
        Array<OneD, Array<OneD, NekDouble> > velfields(m_velocity.num_elements());

        for(int i = 0; i < m_velocity.num_elements(); ++i)
        {
            velfields[i] = m_fields[m_velocity[i]]->UpdatePhys();
        }        
        stdVelocity = GetStdVelocity(velfields);
	
        for(int el = 0; el < n_element; ++el)
        {
            int npoints = m_fields[0]->GetExp(el)->GetTotPoints();
            
            tstep[el] =  CFL[el]/(stdVelocity[el]*cLambda*(ExpOrder[el]-1)*(ExpOrder[el]-1));
        }
	
        NekDouble TimeStep = Vmath::Vmin(n_element,tstep,1);
	
        return TimeStep;
    }
    
    NekDouble IncNavierStokes::v_GetTimeStep(int ExpOrder, NekDouble CFL, NekDouble TimeStability)
    {
        Array<OneD, int> ExpOrderList(m_fields[0]->GetExpSize(),ExpOrder);
        Array<OneD, NekDouble> CFLList(m_fields[0]->GetExpSize(),CFL);
        
        return v_GetTimeStep(ExpOrderList,CFLList,TimeStability);
    }
}
