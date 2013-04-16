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

#include <IncNavierStokesSolver/AdvectionTerms/AdjointAdvection.h>
#include <StdRegions/StdSegExp.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>
#include <MultiRegions/DisContField1D.h>
#include <MultiRegions/DisContField2D.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <MultiRegions/DisContField3DHomogeneous1D.h>
#include <MultiRegions/ContField3DHomogeneous2D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>


namespace Nektar
{
    string AdjointAdvection::className = GetAdvectionTermFactory().RegisterCreatorFunction("Adjoint", AdjointAdvection::create);

    /**
     * Constructor. Creates ...
     *
     * \param 
     * \param
     */

   AdjointAdvection::AdjointAdvection(
            const LibUtilities::SessionReaderSharedPtr&        pSession,
            const SpatialDomains::MeshGraphSharedPtr&          pGraph):
        AdvectionTerm(pSession, pGraph)
	{
	}
	

	void AdjointAdvection::v_InitObject()
	{
	    AdvectionTerm::v_InitObject();
		m_boundaryConditions = MemoryManager<SpatialDomains::BoundaryConditions>
		::AllocateSharedPtr(m_session, m_graph);
		
		//Setting parameters for homogeneous problems
		m_HomoDirec       = 0;
        m_useFFT          = false;
        m_HomogeneousType = eNotHomogeneous;
		m_SingleMode	   =false;
		m_HalfMode		   =false;
		m_MultipleModes    =false;
		
        if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
        {
            std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
            m_spacedim          = 3;
			
            if((HomoStr == "HOMOGENEOUS1D")||(HomoStr == "Homogeneous1D")||
               (HomoStr == "1D")||(HomoStr == "Homo1D"))
            {
                m_HomogeneousType = eHomogeneous1D;
                m_LhomZ           = m_session->GetParameter("LZ");
                m_HomoDirec       = 1;
				
				if(m_session->DefinesSolverInfo("ModeType"))
				{
					m_session->MatchSolverInfo("ModeType","SingleMode",m_SingleMode,false);
					m_session->MatchSolverInfo("ModeType","HalfMode",m_HalfMode,false);
					m_session->MatchSolverInfo("ModeType","MultipleModes",m_MultipleModes,false);
				}
				if(m_session->DefinesSolverInfo("ModeType"))
				{
					if(m_SingleMode)
					{
						m_npointsZ=2;
						
					}
					else if(m_HalfMode)
					{
						m_npointsZ=1;
						
					}
					else if(m_MultipleModes)
					{
						m_npointsZ        = m_session->GetParameter("HomModesZ");
					}
					else
					{
						ASSERTL0(false, "SolverInfo ModeType not valid");	
						
						
					}
				}
				else 
				{
					m_session->LoadParameter("HomModesZ",m_npointsZ);
					
				}
				
            }
			
            if((HomoStr == "HOMOGENEOUS2D")||(HomoStr == "Homogeneous2D")||
               (HomoStr == "2D")||(HomoStr == "Homo2D"))
            {
                m_HomogeneousType = eHomogeneous2D;
                m_session->LoadParameter("HomModesY", m_npointsY);
                m_session->LoadParameter("LY",        m_LhomY);
                m_session->LoadParameter("HomModesZ", m_npointsZ);
                m_session->LoadParameter("LZ",        m_LhomZ);
                m_HomoDirec       = 2;
            }
			
            if((HomoStr == "HOMOGENEOUS3D")||(HomoStr == "Homogeneous3D")||
               (HomoStr == "3D")||(HomoStr == "Homo3D"))
            {
                m_HomogeneousType = eHomogeneous3D;
                m_session->LoadParameter("HomModesX",m_npointsX);
                m_session->LoadParameter("LX",       m_LhomX   );
                m_session->LoadParameter("HomModesY",m_npointsY);
                m_session->LoadParameter("LY",       m_LhomY   );
                m_session->LoadParameter("HomModesZ",m_npointsZ);
                m_session->LoadParameter("LZ",       m_LhomZ   );
                m_HomoDirec       = 3;
            }
			
            if(m_session->DefinesSolverInfo("USEFFT"))
            {
                m_useFFT = true;
            }
        }
        else
        {
            m_npointsZ = 1; // set to default value so can use to identify 2d or 3D (homogeneous) expansions
        }
		
		if(m_session->DefinesSolverInfo("PROJECTION"))
        {
            std::string ProjectStr
			= m_session->GetSolverInfo("PROJECTION");
			
            if((ProjectStr == "Continuous")||(ProjectStr == "Galerkin")||
               (ProjectStr == "CONTINUOUS")||(ProjectStr == "GALERKIN"))
            {
                m_projectionType = MultiRegions::eGalerkin;
            }
            else if(ProjectStr == "DisContinuous")
            {
                m_projectionType = MultiRegions::eDiscontinuous;
            }
            else
            {
                ASSERTL0(false,"PROJECTION value not recognised");
            }
        }
        else
        {
            cerr << "Projection type not specified in SOLVERINFO,"
			"defaulting to continuous Galerkin" << endl;
            m_projectionType = MultiRegions::eGalerkin;
        }
		
		SetUpBaseFields(m_graph);
		ASSERTL0(m_session->DefinesFunction("BaseFlow"),
				 "Base flow must be defined for linearised forms.");
		string file = m_session->GetFunctionFilename("BaseFlow", 0);
		
		
		//Periodic base flows
		if(m_session->DefinesParameter("N_slices"))
		{
            m_session->LoadParameter("N_slices",m_slices);
            if(m_slices>1)
            {
				DFT(file,m_slices);
			}
			else
			{
				ASSERTL0(false,"Number of slices must be a positive number greater than 1");
			}
		}
		//Steady base-flow
		else
		{
			m_slices=1;
			
			//BaseFlow from file
			if (m_session->GetFunctionType("BaseFlow", m_session->GetVariable(0))
				== LibUtilities::eFunctionTypeFile)
			{
				ImportFldBase(file,m_graph,1);
				
			}
			//analytic base flow
			else
			{
				int nq = m_base[0]->GetNpoints();
				Array<OneD,NekDouble> x0(nq);
				Array<OneD,NekDouble> x1(nq);
				Array<OneD,NekDouble> x2(nq);
				
				// get the coordinates (assuming all fields have the same
				// discretisation)
				m_base[0]->GetCoords(x0,x1,x2);
				for(unsigned int i = 0 ; i < m_base.num_elements(); i++)
				{
					LibUtilities::EquationSharedPtr ifunc
					= m_session->GetFunction("BaseFlow", i);
					
					ifunc->Evaluate(x0,x1,x2,m_base[i]->UpdatePhys());
					
					m_base[i]->SetPhysState(true);						
					m_base[i]->FwdTrans_IterPerExp(m_base[i]->GetPhys(),
												   m_base[i]->UpdateCoeffs());
				}
				
			}
			
		}
}
	

    AdjointAdvection::~AdjointAdvection()
    {
    }

    
    void AdjointAdvection::SetUpBaseFields(SpatialDomains::MeshGraphSharedPtr &mesh)
    {
		int nvariables = m_session->GetVariables().size();
        int i;
        m_base = Array<OneD, MultiRegions::ExpListSharedPtr>(nvariables);
        if (m_projectionType == MultiRegions::eGalerkin)
        {
            switch (m_expdim)
            {
				case 1:
				{
					if(m_HomogeneousType == eHomogeneous2D)
					{
						const LibUtilities::PointsKey PkeyY(m_npointsY,LibUtilities::eFourierEvenlySpaced);
						const LibUtilities::BasisKey  BkeyY(LibUtilities::eFourier,m_npointsY,PkeyY);
						const LibUtilities::PointsKey PkeyZ(m_npointsZ,LibUtilities::eFourierEvenlySpaced);
						const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,m_npointsZ,PkeyZ);
						
						for(i = 0 ; i < m_base.num_elements(); i++)
						{
							m_base[i] = MemoryManager<MultiRegions::ContField3DHomogeneous2D>
							::AllocateSharedPtr(m_session,BkeyY,BkeyZ,m_LhomY,m_LhomZ,m_useFFT,m_dealiasing,m_graph,m_session->GetVariable(i));
						}
					}
					
					else {
						
						for(i = 0 ; i < m_base.num_elements(); i++)
						{
							m_base[i] = MemoryManager<MultiRegions::ContField1D>
							::AllocateSharedPtr(m_session,mesh,
												m_session->GetVariable(i));
						}
					}
					
				}
					break;
				case 2:
				{
					if(m_HomogeneousType == eHomogeneous1D)
					{
						if(m_SingleMode)
						{
							const LibUtilities::PointsKey PkeyZ(m_npointsZ,LibUtilities::eFourierSingleModeSpaced);
							const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,m_npointsZ,PkeyZ);
							
							for(i = 0 ; i < m_base.num_elements(); i++)
							{								
								m_base[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
								::AllocateSharedPtr(m_session,BkeyZ,m_LhomZ,m_useFFT,m_dealiasing,m_graph,m_session->GetVariable(i));
								
							} 
						}
						else if(m_HalfMode)
						{
							//1 plane field (half mode expansion)
							const LibUtilities::PointsKey PkeyZ(m_npointsZ,LibUtilities::eFourierSingleModeSpaced);
							const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourierHalfModeRe,m_npointsZ,PkeyZ);
							
							for(i = 0 ; i < m_base.num_elements(); i++)
							{																
								m_base[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
								::AllocateSharedPtr(m_session,BkeyZ,m_LhomZ,m_useFFT,m_dealiasing,m_graph,m_session->GetVariable(i));
							} 
							
						}
						else 
						{
							const LibUtilities::PointsKey PkeyZ(m_npointsZ,LibUtilities::eFourierEvenlySpaced);
							const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,m_npointsZ,PkeyZ);
							
							
							for(i = 0 ; i < m_base.num_elements(); i++)
							{
								m_base[i] = MemoryManager<MultiRegions::ContField3DHomogeneous1D>
								::AllocateSharedPtr(m_session,BkeyZ,m_LhomZ,m_useFFT,m_dealiasing,m_graph,m_session->GetVariable(i));
								m_base[i]->SetWaveSpace(false);
							} 
							
						}
					}
					else
					{
						i = 0;
						MultiRegions::ContField2DSharedPtr firstbase =
                        MemoryManager<MultiRegions::ContField2D>
                        ::AllocateSharedPtr(m_session,mesh,
                                            m_session->GetVariable(i));
						m_base[0]=firstbase;
						
						for(i = 1 ; i < m_base.num_elements(); i++)
						{
							m_base[i] = MemoryManager<MultiRegions::ContField2D>
                            ::AllocateSharedPtr(*firstbase,mesh,
                                                m_session->GetVariable(i));
						}
					}
				}
					break;
				case 3:
				{
					if(m_HomogeneousType == eHomogeneous3D)
					{
						ASSERTL0(false,"3D fully periodic problems not implemented yet");
					}
					else
					{
						MultiRegions::ContField3DSharedPtr firstbase =
						MemoryManager<MultiRegions::ContField3D>
						::AllocateSharedPtr(m_session,mesh,
											m_session->GetVariable(0));
						m_base[0] = firstbase;
						
						for(i = 1 ; i < m_base.num_elements(); i++)
						{
							m_base[i] = MemoryManager<MultiRegions::ContField3D>
							::AllocateSharedPtr(*firstbase,mesh,
												m_session->GetVariable(i));
						}
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
					if(m_HomogeneousType == eHomogeneous2D)
                    {
                        const LibUtilities::PointsKey PkeyY(m_npointsY,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyY(LibUtilities::eFourier,m_npointsY,PkeyY);
                        const LibUtilities::PointsKey PkeyZ(m_npointsZ,LibUtilities::eFourierEvenlySpaced);
                        const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,m_npointsZ,PkeyZ);
						
                        for(i = 0 ; i < m_base.num_elements(); i++)
                        {
                            m_base[i] = MemoryManager<MultiRegions::DisContField3DHomogeneous2D>
							::AllocateSharedPtr(m_session,BkeyY,BkeyZ,m_LhomY,m_LhomZ,m_useFFT,m_dealiasing,m_graph,m_session->GetVariable(i));
                        }
                    }
					else 
					{
						for(i = 0 ; i < m_base.num_elements(); i++)
						{
							m_base[i] = MemoryManager<MultiRegions
                            ::DisContField1D>::AllocateSharedPtr(m_session,mesh,
                                                                 m_session->GetVariable(i));
						}
					}
                    break;
                }
				case 2:
                {
					if(m_HomogeneousType == eHomogeneous1D)
                    {
						
						const LibUtilities::PointsKey PkeyZ(m_npointsZ,LibUtilities::eFourierEvenlySpaced);
						const LibUtilities::BasisKey  BkeyZ(LibUtilities::eFourier,m_npointsZ,PkeyZ);
						
						for(i = 0 ; i < m_base.num_elements(); i++)
						{
							m_base[i] = MemoryManager<MultiRegions::DisContField3DHomogeneous1D>
							::AllocateSharedPtr(m_session,BkeyZ,m_LhomZ,m_useFFT,m_dealiasing,m_graph,m_session->GetVariable(i));
						}
						
						
						
					}
					else
					{
						for(i = 0 ; i < m_base.num_elements(); i++)
						{
							m_base[i] = MemoryManager<MultiRegions
							::DisContField2D>::AllocateSharedPtr(m_session, mesh,
                                                                 m_session->GetVariable(i));
						}
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
            SpatialDomains::MeshGraphSharedPtr pGraph, int cnt)
    {
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> FieldDef;
        std::vector<std::vector<NekDouble> > FieldData;
        int nqtot = m_base[0]->GetTotPoints();
	
        //Get Homogeneous
        LibUtilities::Import(pInfile,FieldDef,FieldData);
		
        int nvar = m_session->GetVariables().size();
        int s;
	
        if(m_session->DefinesSolverInfo("HOMOGENEOUS"))
        {
            std::string HomoStr = m_session->GetSolverInfo("HOMOGENEOUS");
        }
	
        // copy FieldData into m_fields
        for(int j = 0; j < nvar; ++j)
        {
            for(int i = 0; i < FieldDef.size(); ++i)
            {
                if ((m_session->DefinesSolverInfo("HOMOGENEOUS") &&
                    (m_session->GetSolverInfo("HOMOGENEOUS")=="HOMOGENEOUS1D" ||
                     m_session->GetSolverInfo("HOMOGENEOUS")=="1D" ||
                     m_session->GetSolverInfo("HOMOGENEOUS")=="Homo1D")) &&
                    nvar==4)
                {
                    // w-component must be ignored and set to zero.
                    if (j != nvar - 2)
                    {
                        // p component (it is 4th variable of the 3D and corresponds 3nd variable of 2D)
                        s = (j == nvar - 1) ? 2 : j;
			
                        //extraction of the 2D
                        m_base[j]->ExtractDataToCoeffs(
                                            FieldDef[i],
                                            FieldData[i],
                                            FieldDef[i]->m_fields[s],
                                            m_base[j]->UpdateCoeffs());
                    }

                    // Put zero on higher modes
                    int ncplane = (m_base[0]->GetNcoeffs()) / m_npointsZ;

                    if (m_npointsZ > 2)
                    {
                        Vmath::Zero(ncplane*(m_npointsZ-2),
                                    &m_base[j]->UpdateCoeffs()[2*ncplane],1);
                    }
                }
                //2D cases and Homogeneous1D Base Flows
                else
                {
                    bool flag = FieldDef[i]->m_fields[j] ==
                                                    m_session->GetVariable(j);

                    ASSERTL1(flag, (std::string("Order of ") + pInfile
                                    + std::string(" data and that defined in "
                                                  "m_boundaryconditions differs")).c_str());
                    
                    m_base[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                                   FieldDef[i]->m_fields[j],
                                                   m_base[j]->UpdateCoeffs());
                }
            }
            
            if(m_SingleMode || m_HalfMode)
            {
                m_base[j]->SetWaveSpace(true);
		
                m_base[j]->BwdTrans(m_base[j]->GetCoeffs(),
                                    m_base[j]->UpdatePhys());
                
                if(m_SingleMode)
                {
                    //copy the bwd into the second plane for single Mode Analysis
                    int ncplane=(m_base[0]->GetNpoints())/m_npointsZ;
                    Vmath::Vcopy(ncplane,&m_base[j]->GetPhys()[0],1,&m_base[j]->UpdatePhys()[ncplane],1);
                }
            }
            else
            {
                m_base[j]->BwdTrans(m_base[j]->GetCoeffs(),
                                    m_base[j]->UpdatePhys());
                
            }
            
        }
	
        //std::string outname ="BaseFlow.bse";
        //WriteFldBase(outname);
	
        if(m_session->DefinesParameter("N_slices"))
        {
            m_nConvectiveFields = m_base.num_elements()-1;
            
            for(int i=0; i<m_nConvectiveFields;++i)
            {
                
                Vmath::Vcopy(nqtot, &m_base[i]->GetPhys()[0], 1, &m_interp[i][cnt*nqtot], 1);				
            }
            
        }
    }
        
   
    //Evaluation of the advective terms
    void AdjointAdvection::v_ComputeAdvectionTerm(
            Array<OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &pVelocity,
            const Array<OneD, const NekDouble> &pU,
            Array<OneD, NekDouble> &pOutarray,
            int pVelocityComponent,
			NekDouble m_time,
            Array<OneD, NekDouble> &pWk)
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
	
        
		grad0 = Array<OneD, NekDouble> (nPointsTot);
		grad_base_u0 = Array<OneD, NekDouble> (nPointsTot);
		grad_base_v0 = Array<OneD, NekDouble> (nPointsTot);
		grad_base_w0 = Array<OneD, NekDouble> (nPointsTot);	
		
		
		//Evaluation of the base flow for periodic cases
		//(it requires fld files)
		
		if(m_slices>1)
		{				
			if (m_session->GetFunctionType("BaseFlow", 0)
				== LibUtilities::eFunctionTypeFile)
			{
				for(int i=0; i<m_nConvectiveFields;++i)
				{
					UpdateBase(m_slices,m_interp[i],m_base[i]->UpdatePhys(),m_time,m_period);
				}
			}
			else 
			{
				ASSERTL0(false, "Periodic Base flow requires .fld files");	
			}
		}
	
		//Evaluate the Adjoint advection term
        switch(ndim) 
        {
					// 1D
				case 1:
					pFields[0]->PhysDeriv(pU,grad0);
					pFields[0]->PhysDeriv(m_base[0]->GetPhys(),grad_base_u0);
					//Evaluate  U du'/dx
					Vmath::Vmul(nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
					//Evaluate U du'/dx+ u' dU/dx
					Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
					break;
					
					//2D
				case 2:
					
					grad1 = Array<OneD, NekDouble> (nPointsTot);
					grad_base_u1 = Array<OneD, NekDouble> (nPointsTot);
					grad_base_v1 = Array<OneD, NekDouble> (nPointsTot);
				
					pFields[0]->PhysDeriv(pU,grad0,grad1);
				
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
					
					grad1 = Array<OneD, NekDouble> (nPointsTot);
					grad2 = Array<OneD, NekDouble> (nPointsTot);
					grad_base_u1 = Array<OneD, NekDouble> (nPointsTot);
					grad_base_v1 = Array<OneD, NekDouble> (nPointsTot);
					grad_base_w1 = Array<OneD, NekDouble> (nPointsTot);
					
					grad_base_u2 = Array<OneD, NekDouble> (nPointsTot);
					grad_base_v2 = Array<OneD, NekDouble> (nPointsTot);
					grad_base_w2 = Array<OneD, NekDouble> (nPointsTot);
				
					m_base[0]->PhysDeriv(m_base[0]->GetPhys(), grad_base_u0, grad_base_u1,grad_base_u2);
					m_base[0]->PhysDeriv(m_base[1]->GetPhys(), grad_base_v0, grad_base_v1,grad_base_v2);
					m_base[0]->PhysDeriv(m_base[2]->GetPhys(), grad_base_w0, grad_base_w1, grad_base_w2);	
				
					//HalfMode has W(x,y,t)=0
					if(m_HalfMode)
					{
						for(int i=0; i<grad_base_u2.num_elements();++i)
						{
							grad_base_u2[i]=0;
							grad_base_v2[i]=0;
							grad_base_w2[i]=0;
						
						}
					}
				
				pFields[0]->PhysDeriv(pU, grad0, grad1, grad2);
					
				switch (pVelocityComponent)
				{
						//x-equation	
					case 0:
						//Evaluate U du'/dx
						Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
						//Evaluate U du'/dx+ V du'/dy
						Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate U du'/dx+ V du'/dy+W du'/dz
						Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)
						Vmath::Neg(nPointsTot,pOutarray,1);
						//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)+u' dU/dx
						Vmath::Vvtvp(nPointsTot,grad_base_u0,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
						//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)+u'dU/dx+ v' dV/dx
						Vmath::Vvtvp(nPointsTot,grad_base_v0,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
						//Evaluate -(U du'/dx+ V du'/dy+W du'/dz)+u'dU/dx+ v' dV/dx+ w' dW/dz
						Vmath::Vvtvp(nPointsTot,grad_base_w0,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
						break;
						//y-equation	
					case 1:
						//Evaluate U dv'/dx
						Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
						//Evaluate U dv'/dx+ V dv'/dy
						Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate U dv'/dx+ V dv'/dy+W dv'/dz
						Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate -(U dv'/dx+ V dv'/dy+W dv'/dz)
						Vmath::Neg(nPointsTot,pOutarray,1);
						//Evaluate  -(U dv'/dx+ V dv'/dy+W dv'/dz)+u' dU/dy
						Vmath::Vvtvp(nPointsTot,grad_base_u1,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
						//Evaluate  -(U dv'/dx+ V dv'/dy+W dv'/dz)+u' dU/dy +v' dV/dy
						Vmath::Vvtvp(nPointsTot,grad_base_v1,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
						//Evaluate  -(U dv'/dx+ V dv'/dy+W dv'/dz)+u' dU/dy +v' dV/dy+ w' dW/dy
						Vmath::Vvtvp(nPointsTot,grad_base_w1,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
						break;
						
						//z-equation	
					case 2:
						//Evaluate U dw'/dx
						Vmath::Vmul (nPointsTot,grad0,1,m_base[0]->GetPhys(),1,pOutarray,1);
						//Evaluate U dw'/dx+ V dw'/dx
						Vmath::Vvtvp(nPointsTot,grad1,1,m_base[1]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate U dw'/dx+ V dw'/dx+ W dw'/dz
						Vmath::Vvtvp(nPointsTot,grad2,1,m_base[2]->GetPhys(),1,pOutarray,1,pOutarray,1);
						//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)
						Vmath::Neg(nPointsTot,pOutarray,1);
						//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)+u' dU/dz
						Vmath::Vvtvp(nPointsTot,grad_base_u2,1,pVelocity[0],1,pOutarray,1,pOutarray,1);
						//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)+u' dU/dz+v'dV/dz
						Vmath::Vvtvp(nPointsTot,grad_base_v2,1,pVelocity[1],1,pOutarray,1,pOutarray,1);
						//Evaluate -(U dw'/dx+ V dw'/dx+ W dw'/dz)+u' dU/dz+v'dV/dz + w' dW/dz
						Vmath::Vvtvp(nPointsTot,grad_base_w2,1,pVelocity[2],1,pOutarray,1,pOutarray,1);
						break;
				}
					break;
					
            
        default:
            ASSERTL0(false,"dimension unknown");
        }
    }
		
		void AdjointAdvection::UpdateBase( const NekDouble m_slices,
											 Array<OneD, const NekDouble> &inarray,
											 Array<OneD, NekDouble> &outarray,
											 const NekDouble m_time,
											 const NekDouble m_period)
		{
			
			int npoints=m_base[0]->GetTotPoints();
			
			NekDouble BetaT=2*M_PI*fmod (m_time, m_period) / m_period;
			NekDouble phase;
			Array<OneD, NekDouble> auxiliary(npoints);
			
			Vmath::Vcopy(npoints,&inarray[0],1,&outarray[0],1);
			Vmath::Svtvp(npoints, cos(0.5*m_slices*BetaT),&inarray[npoints],1,&outarray[0],1,&outarray[0],1);
			
			for (int i = 2; i < m_slices; i += 2) 
			{
				phase = (i>>1) * BetaT;
				
				Vmath::Svtvp(npoints, cos(phase),&inarray[i*npoints],1,&outarray[0],1,&outarray[0],1);
				Vmath::Svtvp(npoints, sin(phase), &inarray[(i+1)*npoints], 1, &outarray[0], 1,&outarray[0],1);
			}
			
		}
		
		DNekBlkMatSharedPtr AdjointAdvection::GetFloquetBlockMatrix(FloquetMatType mattype, bool UseContCoeffs) const
		{
			DNekMatSharedPtr    loc_mat;
			DNekBlkMatSharedPtr BlkMatrix;
			int n_exp = 0;
			
			n_exp = m_base[0]->GetTotPoints(); // will operatore on m_phys
			
			Array<OneD,unsigned int> nrows(n_exp);
			Array<OneD,unsigned int> ncols(n_exp);
			
			nrows = Array<OneD, unsigned int>(n_exp,m_slices);
			ncols = Array<OneD, unsigned int>(n_exp,m_slices);
			
			MatrixStorage blkmatStorage = eDIAGONAL;
			BlkMatrix = MemoryManager<DNekBlkMat>
			::AllocateSharedPtr(nrows,ncols,blkmatStorage);
			
			
			const LibUtilities::PointsKey Pkey(m_slices,LibUtilities::eFourierEvenlySpaced);
			const LibUtilities::BasisKey  BK(LibUtilities::eFourier,m_slices,Pkey);
			StdRegions::StdSegExp StdSeg(BK);
			
			StdRegions::StdMatrixKey matkey(StdRegions::eFwdTrans,
                                                        StdSeg.DetShapeType(),
                                                        StdSeg);
			
			loc_mat = StdSeg.GetStdMatrix(matkey);
			
			// set up array of block matrices.
			for(int i = 0; i < n_exp; ++i)
			{
				BlkMatrix->SetBlock(i,i,loc_mat);
			}
			
			return BlkMatrix;
		}
		
		//Discrete Fourier Transform for Floquet analysis
		void AdjointAdvection::DFT(const string file, const NekDouble m_slices)
		{
			int npoints=m_base[0]->GetTotPoints();
			
			//Convected fields
			int ConvectedFields=m_base.num_elements()-1;
			
			m_interp= Array<OneD, Array<OneD, NekDouble> > (ConvectedFields);
			for(int i=0; i<ConvectedFields;++i)
			{
				m_interp[i]=Array<OneD,NekDouble>(npoints*m_slices);
			}
			
			//Import the slides into auxiliary vector
			//The base flow should be stored in the form filename_i.bse
			for (int i=0; i< m_slices; ++i)
			{
				char chkout[16] = "";
				sprintf(chkout, "%d", i);
				ImportFldBase(file+"_"+chkout+".bse",m_graph,i);
			} 
			
			
			// Discrete Fourier Transform of the fields
			for(int k=0; k< ConvectedFields;++k)
			{
#ifdef NEKTAR_USING_FFTW
				
				//Discrete Fourier Transform using FFTW
				
				
				Array<OneD, NekDouble> fft_in(npoints*m_slices);
				Array<OneD, NekDouble> fft_out(npoints*m_slices);
				
				Array<OneD, NekDouble> m_tmpIN(m_slices);
				Array<OneD, NekDouble> m_tmpOUT(m_slices);
				
				//Shuffle the data
				for(int j= 0; j < m_slices; ++j)
				{
					Vmath::Vcopy(npoints,&m_interp[k][j*npoints],1,&(fft_in[j]),m_slices);
				}
				
				m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", m_slices);
				
				//FFT Transform
				for(int i=0; i<npoints; i++)
				{
					m_FFT->FFTFwdTrans(m_tmpIN =fft_in + i*m_slices, m_tmpOUT =fft_out + i*m_slices);
					
				}
				
				//Reshuffle data
				for(int s = 0; s < m_slices; ++s)
				{						
					Vmath::Vcopy(npoints,&fft_out[s],m_slices,&m_interp[k][s*npoints],1);
					
				}
				
				Vmath::Zero(fft_in.num_elements(),&fft_in[0],1);
				Vmath::Zero(fft_out.num_elements(),&fft_out[0],1);				
#else
				//Discrete Fourier Transform using MVM
				
				
				DNekBlkMatSharedPtr blkmat;
				blkmat = GetFloquetBlockMatrix(eForwardsPhys);
				
				int nrows = blkmat->GetRows();
				int ncols = blkmat->GetColumns();
				
				Array<OneD, NekDouble> sortedinarray(ncols);
				Array<OneD, NekDouble> sortedoutarray(nrows);
				
				//Shuffle the data
				for(int j= 0; j < m_slices; ++j)
				{
					Vmath::Vcopy(npoints,&m_interp[k][j*npoints],1,&(sortedinarray[j]),m_slices);
				}
				
				// Create NekVectors from the given data arrays
				NekVector<NekDouble> in (ncols,sortedinarray,eWrapper);
				NekVector<NekDouble> out(nrows,sortedoutarray,eWrapper);
				
				// Perform matrix-vector multiply.
				out = (*blkmat)*in;
				
				//Reshuffle data
				for(int s = 0; s < m_slices; ++s)
				{	
					Vmath::Vcopy(npoints,&sortedoutarray[s],m_slices,&m_interp[k][s*npoints],1);
				}
                                
				Vmath::Zero(sortedinarray.num_elements(),&sortedinarray[0],1);
				Vmath::Zero(sortedoutarray.num_elements(),&sortedoutarray[0],1);	
								
#endif
				
				//scaling of the Fourier coefficients
				NekDouble j=-1;
				for (int i = 2; i < m_slices; i += 2) 
				{
					Vmath::Smul(2*npoints,j,&m_interp[k][i*npoints],1,&m_interp[k][i*npoints],1);
					j=-j;
					
				}
				
			}
			
			if(m_session->DefinesParameter("period"))
			{
				m_period=m_session->GetParameter("period");
			}
			else 
			{
				m_period=(m_session->GetParameter("TimeStep")*m_slices)/(m_slices-1.);
			}	
			
			
		}
		
	
} //end of namespace

