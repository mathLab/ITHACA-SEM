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
#include <cstdio>
#include <cstdlib>

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
		SetUpBaseFields(m_graph);
		ASSERTL0(m_session->DefinesFunction("BaseFlow"),
				 "Base flow must be defined for linearised forms.");
		ASSERTL0(m_session->GetFunctionType("BaseFlow")
				 == LibUtilities::eFunctionTypeFile,
				 "Base flow must be provided in a file.");
		string file = m_session->GetFunctionFilename("BaseFlow");
		
		//Periodic base flows
		if(m_session->DefinesParameter("N_slices"))
		{
		    m_session->LoadParameter("N_slices", m_slices);

			if(m_slices>1)
			{
				
				int npoints=m_base[0]->GetTotPoints();
				Array<OneD, NekDouble> fft_in(npoints*m_slices);
				Array<OneD, NekDouble> fft_out(npoints*m_slices);
				
				Array<OneD, NekDouble> m_tmpIN(m_slices);
				Array<OneD, NekDouble> m_tmpOUT(m_slices);
				
				//Convected fields
				int ConvectedFields=m_base.num_elements()-1;
				
				m_interp= Array<OneD, Array<OneD, NekDouble> > (ConvectedFields);
				for(int i=0; i<ConvectedFields;++i)
				{
					m_interp[i]=Array<OneD,NekDouble>(npoints*m_slices);
				}
				
				cout << "file " << file << endl;
				
				//Import the slides into auxiliary vector
				for (int i=0; i< m_slices; ++i)
				{
					char chkout[16] = "";
					sprintf(chkout, "%d", i);
					ImportFldBase(file+"_"+chkout+".fld",m_graph,i);
				} 
				
				m_useFFTW=false;
				if(m_session->DefinesSolverInfo("USEFFT"))
				{
					m_useFFTW = true;
				}
				
				//Factory for FFT transformation
				if(m_useFFTW)
				{
					m_FFT = LibUtilities::GetNektarFFTFactory().CreateInstance("NekFFTW", m_slices);
				}
				else 
				{
					ASSERTL0(false, "Time interpolation not implemented");
				}
				
				// Discrete Fourier Transform of the fields
				for(int k=0; k< ConvectedFields;++k)
				{
					//Shuffle the data
					for(int j= 0; j < m_slices; ++j)
					{
						Vmath::Vcopy(npoints,&m_interp[k][j*npoints],1,&(fft_in[j]),m_slices);
					}
					
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
					
					for(int r=0; r<fft_in.num_elements();++r)
					{
						fft_in[0]=0;
						fft_out[0]=0;
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
			else{
				
				ASSERTL0(false,"Number of slices must be a positive number");
			}
		}
		//Steady base-flow
		else
		{
			m_slices=1;
			ImportFldBase(file,m_graph,1);
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
        
        if (m_projectionType = MultiRegions::eGalerkin)
        {
            switch (m_expdim)
            {
            case 1:
	        {
                    for(i = 0 ; i < m_base.num_elements(); i++)
                    {
                        m_base[i] = MemoryManager<MultiRegions::ContField1D>
                            ::AllocateSharedPtr(m_session,mesh,
                                                m_session->GetVariable(i));
                    }
	        }
                break;
            case 2:
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
                break;
            case 3:
	        {
                    MultiRegions::ContField3DSharedPtr firstbase =
                        MemoryManager<MultiRegions::ContField3D>
                        ::AllocateSharedPtr(m_session,mesh,
                                            m_session->GetVariable(i));
                    m_base[0] = firstbase;
                    
                    for(i = 1 ; i < m_base.num_elements(); i++)
                    {
                        m_base[i] = MemoryManager<MultiRegions::ContField3D>
                            ::AllocateSharedPtr(*firstbase,mesh,
                                                m_session->GetVariable(i));
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
                    for(i = 0 ; i < m_base.num_elements(); i++)
                    {
                        m_base[i] = MemoryManager<MultiRegions
                            ::DisContField1D>::AllocateSharedPtr(m_session,mesh,
                                                                 m_session->GetVariable(i));
                    }
                    break;
                }
            case 2:
                {
                    for(i = 0 ; i < m_base.num_elements(); i++)
                    {
                        m_base[i] = MemoryManager<MultiRegions
                            ::DisContField2D>::AllocateSharedPtr(m_session, mesh,
                                                                 m_session->GetVariable(i));
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
        std::vector<SpatialDomains::FieldDefinitionsSharedPtr> FieldDef;
        std::vector<std::vector<NekDouble> > FieldData;
		int numfields=m_base.num_elements();
		int nqtot = m_base[0]->GetTotPoints();

        pGraph->Import(pInfile,FieldDef,FieldData);

        int nvar = m_session->GetVariables().size();

        // copy FieldData into m_fields
        for(int j = 0; j < nvar; ++j)
        {
            for(int i = 0; i < FieldDef.size(); ++i)
            {
                bool flag = FieldDef[i]->m_fields[j]
                    == m_session->GetVariable(j);
                ASSERTL1(flag, (std::string("Order of ") + pInfile
                                + std::string(" data and that defined in "
                                              "m_boundaryconditions differs")).c_str());
                
                m_base[j]->ExtractDataToCoeffs(FieldDef[i], FieldData[i],
                                               FieldDef[i]->m_fields[j]);
            }
            m_base[j]->BwdTrans(m_base[j]->GetCoeffs(),
                                m_base[j]->UpdatePhys());
		}
			
			if(m_session->DefinesParameter("N_slices"))
			{
				
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
	
        
        grad0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_u0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_v0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
        grad_base_w0 = Array<OneD, NekDouble> (ndim*nPointsTot, 0.0);
		
		//Evaluation of the base flow for periodic cases
		if(m_slices>1)
		{
			for(int i=0; i<m_nConvectiveFields;++i)
			{
				UpdateBase(m_slices,m_interp[i],m_base[i]->UpdatePhys(),m_time,m_period);
			}
		}
	
		//Evaluate the Adjoint advection term
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
				Vmath::Svtvp(npoints, -sin(phase), &inarray[(i+1)*npoints], 1, &outarray[0], 1,&outarray[0],1);
			}
			
		}
		
	
} //end of namespace

