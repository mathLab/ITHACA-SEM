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

#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>

namespace Nektar
{
	/*
    NekDouble IncNavierStokes::v_GetTimeStep(const Array<OneD,int> ExpOrder, const Array<OneD,NekDouble> CFL, NekDouble timeCFL)
    { 
        
        int nvariables      = m_fields.num_elements();
        int nTotQuadPoints  = GetTotPoints();
        int n_element       = m_fields[0]->GetExpSize(); 
        
        const NekDouble minLengthStdTri  = 1.414213;
        const NekDouble minLengthStdQuad = 2.0;
        const NekDouble cLambda          = 0.2; // Spencer book pag. 317
        
        Array<OneD, NekDouble> tstep      (n_element, 0.0);
        Array<OneD, NekDouble> stdVelocity(n_element, 0.0);
        Array<OneD, Array<OneD, NekDouble> > velfields(
                                                m_velocity.num_elements());

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
    
    Array<OneD,NekDouble> IncNavierStokes::GetStdVelocity(const Array<OneD, Array<OneD,NekDouble> > inarray)
	{
        // Checking if the problem is 2D
        ASSERTL0(m_expdim>=2,"Method not implemented for 1D");
        
        int nTotQuadPoints  = GetTotPoints();
        int n_element  = m_fields[0]->GetExpSize();       // number of element in the mesh
        int nvel = inarray.num_elements();
        int npts = 0;
        
        NekDouble pntVelocity;
        
        // Getting the standard velocity vector on the 2D normal space
        Array<OneD, Array<OneD, NekDouble> > stdVelocity(nvel);
		
        Array<OneD, NekDouble> stdV(n_element,0.0);
        for (int i = 0; i < nvel; ++i)
        {
            stdVelocity[i] = Array<OneD, NekDouble>(nTotQuadPoints);
        }
		
        if(nvel == 2)
        {
            for(int el = 0; el < n_element; ++el)
            { 
                
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                
                Array<OneD, const NekDouble> jac  = m_fields[0]->GetExp(el)->GetGeom2D()->GetJac();
                Array<TwoD, const NekDouble> gmat = m_fields[0]->GetExp(el)->GetGeom2D()->GetGmat();
                
                if(m_fields[0]->GetExp(el)->GetGeom2D()->GetGtype() == SpatialDomains::eDeformed)
                {
                    for(int i=0; i<n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i] + gmat[2][i]*inarray[1][i];
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i] + gmat[3][i]*inarray[1][i];
                    }
                }
                else
                {
                    for(int i=0; i<n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i] + gmat[2][0]*inarray[1][i];
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i] + gmat[3][0]*inarray[1][i];
                    }
                }
                
                
                for(int i=0; i<n_points; i++)
                {
                    pntVelocity = sqrt(stdVelocity[0][i]*stdVelocity[0][i] + stdVelocity[1][i]*stdVelocity[1][i]);
                    if(pntVelocity>stdV[el])
                    {
                        stdV[el] = pntVelocity;
                    }
                    
                }
            }
        }
        else
        {
            for(int el = 0; el < n_element; ++el)
            { 
                
                int n_points = m_fields[0]->GetExp(el)->GetTotPoints();
                
                Array<OneD, const NekDouble> jac  = m_fields[0]->GetExp(el)->GetGeom3D()->GetJac();
                Array<TwoD, const NekDouble> gmat = m_fields[0]->GetExp(el)->GetGeom3D()->GetGmat();
                
                if(m_fields[0]->GetExp(el)->GetGeom3D()->GetGtype() == SpatialDomains::eDeformed)
                {
                    for(int i=0; i<n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][i]*inarray[0][i] + gmat[3][i]*inarray[1][i] + gmat[6][i]*inarray[2][i];
                        stdVelocity[1][i] = gmat[1][i]*inarray[0][i] + gmat[4][i]*inarray[1][i] + gmat[7][i]*inarray[2][i];
                        stdVelocity[2][i] = gmat[2][i]*inarray[0][i] + gmat[5][i]*inarray[1][i] + gmat[8][i]*inarray[2][i];
                    }
                }
                else
                {
                    Array<OneD, const NekDouble> jac  = m_fields[0]->GetExp(el)->GetGeom3D()->GetJac();
                    Array<TwoD, const NekDouble> gmat = m_fields[0]->GetExp(el)->GetGeom3D()->GetGmat();
                    
                    for(int i=0; i<n_points; i++)
                    {
                        stdVelocity[0][i] = gmat[0][0]*inarray[0][i] + gmat[3][0]*inarray[1][i] + gmat[6][0]*inarray[2][i];
                        stdVelocity[1][i] = gmat[1][0]*inarray[0][i] + gmat[4][0]*inarray[1][i] + gmat[7][0]*inarray[2][i];
                        stdVelocity[2][i] = gmat[2][0]*inarray[0][i] + gmat[5][0]*inarray[1][i] + gmat[8][0]*inarray[2][i];
                    }
                }
                
                for(int i=0; i<n_points; i++)
                {
                    pntVelocity = sqrt(stdVelocity[0][i]*stdVelocity[0][i] + stdVelocity[1][i]*stdVelocity[1][i] + stdVelocity[2][i]*stdVelocity[2][i]);
                    if(pntVelocity>stdV[el])
                    {
                        stdV[el] = pntVelocity;
                    }
                }
            }
        }
		
        return stdV;
	}
     */
}
