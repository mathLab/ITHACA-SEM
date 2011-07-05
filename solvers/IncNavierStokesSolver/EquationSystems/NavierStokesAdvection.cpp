///////////////////////////////////////////////////////////////////////////////
//
// File NavierStokesAdvection.cpp
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
// Description: Evaluation of the Navier Stokes advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/NavierStokesAdvection.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
    /**
     * Basic construnctor
     */
    NavierStokesAdvection::NavierStokesAdvection(void):
        AdvectionTerm()
    {     
    }
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

    NavierStokesAdvection::NavierStokesAdvection(
            LibUtilities::CommSharedPtr                 pComm,
            LibUtilities::SessionReaderSharedPtr        pSession,
            SpatialDomains::MeshGraphSharedPtr          pGraph,
            SpatialDomains::BoundaryConditionsSharedPtr pBoundaryConditions):
        AdvectionTerm(pComm, pSession, pGraph, pBoundaryConditions)
	
    {

	}
	
	NavierStokesAdvection::~NavierStokesAdvection()
	{
	}
	
	//Advection function
	void NavierStokesAdvection:: v_DoAdvection(	Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
											    const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
												Array<OneD, Array<OneD, NekDouble> > &pOutarray,
											    Array<OneD, NekDouble> &pWk)
	
	{
		int i,j;
		int VelDim;
		int numfields = pFields.num_elements();
		std::string velids[] = {"u","v","w"};
        int nqtot      = pFields[0]->GetTotPoints();

	    // Assume all fields but last to be convected by velocity. 
		m_nConvectiveFields=numfields-1;
        
		m_velocity = Array<OneD,int>(m_nConvectiveFields);
		
		for(i = 0; i <m_nConvectiveFields; ++i)
        {
            for(j = 0; j < numfields; ++j)
            {
                std::string var = m_boundaryConditions->GetVariable(j);
                if(NoCaseStringCompare(velids[i],var) == 0)
                {
                    m_velocity[i] = j;
                    break;
                }
                
                if(j == numfields)
                {
                    std::string error = "Failed to find field: " + var; 
                    ASSERTL0(false,error.c_str());
                }
            }
        }
		
        VelDim     = m_velocity.num_elements();
        
		Array<OneD, Array<OneD, NekDouble> > velocity(VelDim);
        
		Array<OneD, NekDouble > Deriv;
        
        for(i = 0; i < VelDim; ++i)
        {
            velocity[i] = pInarray[m_velocity[i]];
        }
		
		// Set up Derivative work space; 
        if(pWk.num_elements())
        {
            ASSERTL0(pWk.num_elements() > nqtot*VelDim,"Workspace is not sufficient");
            Deriv = pWk;
        }
        else
        {
            Deriv = Array<OneD, NekDouble> (nqtot*VelDim);
        }
		
		for(i=0; i< m_nConvectiveFields; ++i)
		{
		 ComputeAdvectionTerm(m_boundaryConditions, pFields,velocity,pInarray[i],pOutarray[i],Deriv);
		 
			
/**			if(i == 0)
			{
			   for (int k=0; k<nqtot;++k)
			   {
				   pOutarray[0][k] = -2;
			   }
			}
*/			
		Vmath::Neg(nqtot,pOutarray[i],1);
			
			
		}
	 }


	//Evaluation of the advective terms
    void NavierStokesAdvection::ComputeAdvectionTerm(SpatialDomains::BoundaryConditionsSharedPtr &pBoundaryConditions,
											Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
											const Array<OneD, Array<OneD, NekDouble> > &pV,
											const Array<OneD, const NekDouble> &pU,
										     Array<OneD, NekDouble> &pOutarray,
											 Array<OneD, NekDouble> &pWk)
    {
		// use dimension of Velocity vector to dictate dimension of operation
        int ndim       = pV.num_elements();
		
        // ToDo: here we should add a check that V has right dimension
		
        int nPointsTot = pFields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2;
		
        grad0 = Array<OneD, NekDouble> (nPointsTot);
		
        // Evaluate V\cdot Grad(u)
        switch(ndim)
        {
			case 1:
				pFields[0]->PhysDeriv(pU,grad0);
				Vmath::Vmul(nPointsTot,grad0,1,pV[0],1,pOutarray,1);
				break;
			case 2:
				grad1 = Array<OneD, NekDouble> (nPointsTot);
		        pFields[0]->PhysDeriv(pU,grad0,grad1);
				Vmath::Vmul (nPointsTot,grad0,1,pV[0],1,pOutarray,1);
			    Vmath::Vvtvp(nPointsTot,grad1,1,pV[1],1,pOutarray,1,pOutarray,1);
				break;	 
			case 3:
				grad1 = Array<OneD, NekDouble> (nPointsTot);
				grad2 = Array<OneD, NekDouble> (nPointsTot);
				pFields[0]->PhysDeriv(pU,grad0,grad1,grad2);
				Vmath::Vmul(nPointsTot,grad0,1,pV[0],1,pOutarray,1);
				Vmath::Vvtvp(nPointsTot,grad1,1,pV[1],1,pOutarray,1,pOutarray,1);
				Vmath::Vvtvp(nPointsTot,grad2,1,pV[2],1,pOutarray,1,pOutarray,1);
				break;
			default:
				ASSERTL0(false,"dimension unknown");
        }
	
			
	}

} //end of namespace

