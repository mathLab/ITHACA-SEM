///////////////////////////////////////////////////////////////////////////////
//
// File AdjointAdvection.cpp
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
// Description: Evaluation of the adjoint advective term
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/AdjointAdvection.h>
#include <cstdio>
#include <cstdlib>

namespace Nektar
{
    /**
     * Basic construnctor
     */
    AdjointAdvection::AdjointAdvection(void):
        AdvectionTerm()
    {     
    }
    
    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

    AdjointAdvection::AdjointAdvection(
            LibUtilities::CommSharedPtr                 pComm,
            LibUtilities::SessionReaderSharedPtr        pSession,
            SpatialDomains::MeshGraphSharedPtr          pGraph,
            SpatialDomains::BoundaryConditionsSharedPtr pBoundaryConditions):
        AdvectionTerm(pComm, pSession, pGraph, pBoundaryConditions)
	{
        SetUpBaseFields(pGraph);
		   ImportFldBase(pSession->GetFilename().substr(0,pSession->GetFilename().find_last_of('.')) + ".bse",pGraph,pBoundaryConditions);
	}
	

	AdjointAdvection::~AdjointAdvection()
	{
	}
	

	void AdjointAdvection::SetUpBaseFields(SpatialDomains::MeshGraphSharedPtr &mesh)
	{
	    int nvariables = m_boundaryConditions->GetNumVariables();
	    int i;
	    m_base = Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);

	    if (m_projectionType = ADRBase::eGalerkin)
	    {
	        switch (m_expdim)
	        {
	        case 1:
	        {
                SpatialDomains::MeshGraph1DSharedPtr mesh1D;

                if( !(mesh1D = boost::dynamic_pointer_cast<
                      SpatialDomains::MeshGraph1D>(mesh)) )
                {
                    ASSERTL0(false,"Dynamics cast failed");
                }

                for(i = 0 ; i < m_base.num_elements(); i++)
                {
                    m_base[i] = MemoryManager<MultiRegions::ContField1D>
                                ::AllocateSharedPtr(m_comm,*mesh1D,
                                                    *m_boundaryConditions,i);
                }
	        }
            break;
	        case 2:
	        {
                SpatialDomains::MeshGraph2DSharedPtr mesh2D;

                if(!(mesh2D = boost::dynamic_pointer_cast<
                     SpatialDomains::MeshGraph2D>(mesh)))
                {
                    ASSERTL0(false,"Dynamics cast failed");
                }

                i = 0;
                MultiRegions::ContField2DSharedPtr firstbase =
                        MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(m_comm,*mesh2D,
                                                    *m_boundaryConditions,i);
                m_base[0]=firstbase;

                for(i = 1 ; i < m_base.num_elements(); i++)
                {
                    m_base[i] = MemoryManager<MultiRegions::ContField2D>
                                ::AllocateSharedPtr(*firstbase,*mesh2D,
                                                    *m_boundaryConditions,i);
                }
	        }
            break;
	        case 3:
	        {
	            SpatialDomains::MeshGraph3DSharedPtr mesh3D;

                if(!(mesh3D = boost::dynamic_pointer_cast<
                     SpatialDomains::MeshGraph3D>(mesh)))
                {
                    ASSERTL0(false,"Dynamics cast failed");
                }

                MultiRegions::ContField3DSharedPtr firstbase =
                        MemoryManager<MultiRegions::ContField3D>
                                ::AllocateSharedPtr(m_comm,*mesh3D,
                                                    *m_boundaryConditions,i);
                m_base[0] = firstbase;

                for(i = 1 ; i < m_base.num_elements(); i++)
                {
                    m_base[i] = MemoryManager<MultiRegions::ContField3D>
                    ::AllocateSharedPtr(*firstbase,
                                        *mesh3D,*m_boundaryConditions,i);
                }
	        }
            break;
            default:
                ASSERTL0(false,"Expansion dimension not recognised");
                break;
	        }
	    }
	    else
	    {
            switch(m_expdim)
            {
                case 1:
                {
                    SpatialDomains::MeshGraph1DSharedPtr mesh1D;

                    if(!(mesh1D = boost::dynamic_pointer_cast<SpatialDomains
                         ::MeshGraph1D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    for(i = 0 ; i < m_base.num_elements(); i++)
                    {
                        m_base[i] = MemoryManager<MultiRegions
                        ::DisContField1D>::AllocateSharedPtr(m_comm,*mesh1D,
                                                             *m_boundaryConditions,i);
                    }
                    break;
                }
                case 2:
                {
                    SpatialDomains::MeshGraph2DSharedPtr mesh2D;

                    if(!(mesh2D = boost::dynamic_pointer_cast<SpatialDomains
                         ::MeshGraph2D>(mesh)))
                    {
                        ASSERTL0(false,"Dynamics cast failed");
                    }

                    for(i = 0 ; i < m_base.num_elements(); i++)
                    {
                        m_base[i] = MemoryManager<MultiRegions
                        ::DisContField2D>::AllocateSharedPtr(m_comm, *mesh2D,
                                                             *m_boundaryConditions,i);
                    }
                    break;
                }
                case 3:
                    ASSERTL0(false,"3 D not set up");
                default:
                    ASSERTL0(false,"Expansion dimension not recognised");
                    break;
            }
	    }

	}


    /**
     * Import field from infile and load into \a m_fields. This routine will
     * also perform a \a BwdTrans to ensure data is in both the physical and
     * coefficient storage.
     * @param   infile          Filename to read.
     */
    void AdjointAdvection::ImportFldBase(std::string pInfile,
            SpatialDomains::MeshGraphSharedPtr pGraph,
            SpatialDomains::BoundaryConditionsSharedPtr &pBoundaryConditions)
    {
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef;
        std::vector<std::vector<NekDouble> > FieldData;

        pGraph->Import(pInfile,FieldDef,FieldData);

        int nvar = pBoundaryConditions->GetNumVariables();

        // copy FieldData into m_fields
        for(int j = 0; j < nvar; ++j)
        {
            for(int i = 0; i < FieldDef.size(); ++i)
            {
                bool flag = FieldDef[i]->m_fields[j]
                                    == pBoundaryConditions->GetVariable(j);
                ASSERTL1(flag, (std::string("Order of ") + pInfile
                            + std::string(" data and that defined in "
                                    "m_boundaryconditions differs")).c_str());

                m_base[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                                 FieldDef[i]->m_fields[j]);
            }
            m_base[j]->BwdTrans(m_base[j]->GetCoeffs(),
                                  m_base[j]->UpdatePhys());
        }
    }


	void AdjointAdvection:: v_DoAdvection(
										  Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
										  const Array<OneD, const Array<OneD, NekDouble> > &pInarray,
										  Array<OneD, Array<OneD, NekDouble> > &pOutarray,
										  Array<OneD, NekDouble> &pWk)
	
	{
		
		int i,j;
		int VelDim;
        int numfields = pFields.num_elements();
		std::string velids[] = {"u","v","w"};

        int nqtot      = pFields[0]->GetTotPoints();
		 
		nvariables=m_boundaryConditions->GetNumVariables();
		// Assume all fields but last to be convected by velocity. 
		m_nConvectiveFields=numfields-1;

		m_velocity = Array<OneD,int>(m_nConvectiveFields);		
		
		for(i = 0; i < m_expdim; ++i)
        {
            for(j = 0; j < m_nConvectiveFields; ++j)
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
		 
		// Name of of the file where the base flow is read
		std::string filefld =  m_sessionName + "-Base.fld";
    
		for(i = 0; i < m_nConvectiveFields; ++i)
		{
		    //cout <<" ----i ="<< i <<endl;
			ComputeAdvectionTerm(pFields, i, velocity,pOutarray[i]);
		    Vmath::Neg(nqtot,pOutarray[i],1);
		}
	 }


	//Evaluation of the advective terms
    void AdjointAdvection::ComputeAdvectionTerm(
										  Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
										  int pVelocityComponent,
										  const Array<OneD, Array<OneD, NekDouble> > &pVelocity,
										  Array<OneD, NekDouble> &pOutarray
										  )
    {
        int ndim       = m_nConvectiveFields;
        int nPointsTot = pFields[0]->GetNpoints();
        Array<OneD, NekDouble> grad0,grad1,grad2;
		
		//Evaluation of the gradiend of each component of the base flow
		//\nabla U
		Array<OneD, NekDouble> grad_base_u0,grad_base_u1,grad_base_u2;
		// \nabla V
		Array<OneD, NekDouble> grad_base_v0,grad_base_v1,grad_base_v2;
	    // \nabla W
		Array<OneD, NekDouble> grad_base_w0,grad_base_w1,grad_base_w2;
		

        grad0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_u0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_v0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_w0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
		
		//Evaluate the adjoint advection term
		switch(ndim) 
        {
			// 1D
			case 1:
				pFields[0]->PhysDeriv(pVelocity[pVelocityComponent],grad0);
				pFields[0]->PhysDeriv(m_base[0]->GetPhys(),grad_base_u0);
				//Evaluate  U du'/dx
				Vmath::Vmul(nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
				//Evaluate U du'/dx+ u' dU/dx
				Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
				break;
			
			//2D
			case 2:
				
				grad1 = grad0 + nPointsTot;
				grad_base_u1 = grad_base_u0 + nPointsTot;
				grad_base_v1 = grad_base_v0 +nPointsTot;
				
				pFields[0]->PhysDeriv(pVelocity[pVelocityComponent],grad0,grad1);

				//Derivates of the base flow
				pFields[0]-> PhysDeriv(m_base[0]->GetPhys(), grad_base_u0, grad_base_u1);
				pFields[0]-> PhysDeriv(m_base[1]->GetPhys(), grad_base_v0, grad_base_v1);

				//Since the components of the velocity are passed one by one, it is necessary to distinguish which
				//term is consider
				switch (pVelocityComponent)
				{
					//x-equation
                    case 0:
						// Evaluate U du'/dx
                        Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
                        //Evaluate U du'/dx+ V du'/dy
                        Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate - (U du'/dx+ V du'/dy)
						Vmath::Neg(nPointsTot,pOutarray,1);
                        //Evaluate -(U du'/dx+ V du'/dy)+u' dU/dx
                        Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
                        //Evaluate -(U du'/dx+ V du'/dy) +u' dU/dx +v' dV/dx
                        Vmath::Vvtvp(nPointsTot,grad_base_v0,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
						break;
			
                    //y-equation
                    case 1:
                        // Evaluate U dv'/dx
                        Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
                        //Evaluate U dv'/dx+ V dv'/dy
                        Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate -(U dv'/dx+ V dv'/dy)
						Vmath::Neg(nPointsTot,pOutarray,1);
                        //Evaluate (U dv'/dx+ V dv'/dy)+u' dU/dy
                        Vmath::Vvtvp(nPointsTot,grad_base_u1,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
                        //Evaluate (U dv'/dx+ V dv'/dy +u' dv/dx)+v' dV/dy
                        Vmath::Vvtvp(nPointsTot,grad_base_v1,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
						break;
				}
				break;
				
								
			//3D
			case 3:
				
				//NOT YET IMPLEMENTED
				grad1 = grad0 + nPointsTot;
				grad2 = grad1 + nPointsTot;
				grad_base_u1 = grad_base_u0 + nPointsTot;
				grad_base_v1 = grad_base_v0 +nPointsTot;
				grad_base_u2 = grad_base_u1 +nPointsTot;
				grad_base_v2 = grad_base_v1 +nPointsTot;
				
				pFields[0]->PhysDeriv(pVelocity[pVelocityComponent], grad0, grad1, grad2);
				
			    pFields[0]->PhysDeriv(m_base[0]->GetPhys(), grad_base_u0, grad_base_u1,grad_base_u2);
				pFields[0]->PhysDeriv(m_base[1]->GetPhys(), grad_base_v0, grad_base_v1,grad_base_v2);
				pFields[0]->PhysDeriv(m_base[2]->GetPhys(), grad_base_w0, grad_base_w1, grad_base_w2);
				
				switch (pVelocityComponent)
				{
					//x-equation	
                    case 0:
                        //Evaluate U du'/dx
                        Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
                        //Evaluate U du'/dx+ V du'/dy
                        Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
                        //Evaluate (U du'/dx+ V du'/dy)+u' dU/dx
                        Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
                        //Evaluate (U du'/dx+ V du'/dy +u' dU/dx)+v' dU/dy
                        Vmath::Vvtvp(nPointsTot,grad_base_u1,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
                        //Evaluate (U du'/dx+ V du'/dy +u' dU/dx +v' dU/dy) + W du'/dz
                        Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
                        //Evaluate (U du'/dx+ V du'/dy +u' dU/dx +v' dU/dy + W du'/dz)+ w' dU/dz
                        Vmath::Vvtvp(nPointsTot,grad_base_u2,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
						break;
					//y-equation	
                    case 1:
                        //Evaluate U dv'/dx
                        Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
                        //Evaluate U dv'/dx+ V dv'/dy
                        Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
                        //Evaluate (U dv'/dx+ V dv'/dy)+u' dV/dx
                        Vmath::Vvtvp(nPointsTot,grad_base_v0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
                        //Evaluate (U du'/dx+ V du'/dy +u' dV/dx)+v' dV/dy
                        Vmath::Vvtvp(nPointsTot,grad_base_v1,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
                        //Evaluate (U du'/dx+ V dv'/dy +u' dV/dx +v' dV/dy) + W du'/dz
                        Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
                        //Evaluate (U du'/dx+ V dv'/dy +u' dV/dx +v' dV/dy + W dv'/dz)+ w' dV/dz
                        Vmath::Vvtvp(nPointsTot,grad_base_v2,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
						break;
 
					//z-equation	
                    case 2:
                        //Evaluate U dw'/dx
                        Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
                        //Evaluate U dw'/dx+ V dw'/dx
                        Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
                        //Evaluate (U dw'/dx+ V dw'/dx)+u' dW/dx
                        Vmath::Vvtvp(nPointsTot,grad_base_w0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
                        //Evaluate (U dw'/dx+ V dw'/dx +w' dW/dx)+v' dW/dy
                        Vmath::Vvtvp(nPointsTot,grad_base_w1,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
                        //Evaluate (U dw'/dx+ V dw'/dx +u' dW/dx +v' dW/dy) + W dw'/dz
                        Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
                        //Evaluate (U dw'/dx+ V dw'/dx +u' dW/dx +v' dW/dy + W dw'/dz)+ w' dW/dz
                        Vmath::Vvtvp(nPointsTot,grad_base_w2,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
						break;
                }
				break;
				
            default:
                ASSERTL0(false,"dimension unknown");
        }
    }
	
	
} //end of namespace

