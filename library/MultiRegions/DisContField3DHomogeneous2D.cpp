//////////////////////////////////////////////////////////////////////////////
//
// File DisContField3DHomogeneous2D.cpp
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
// Description: Field definition for 3D domain with boundary
// conditions using LDG flux and a 2D homogeneous directions
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/ExpList1DHomogeneous2D.h>
#include <MultiRegions/DisContField3DHomogeneous2D.h>

namespace Nektar
{
    namespace MultiRegions
    {

        DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(void):
            ExpList3DHomogeneous2D(),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(LibUtilities::SessionReaderSharedPtr &pSession,
                                                                 const LibUtilities::BasisKey &HomoBasis_y,
																 const LibUtilities::BasisKey &HomoBasis_z,
                                                                 const NekDouble lhom_y,
																 const NekDouble lhom_z,
																 bool useFFT):
            ExpList3DHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(const DisContField3DHomogeneous2D &In, bool DeclareLinesSetCoeffPhys):
            ExpList3DHomogeneous2D (In,false),
            m_bndCondExpansions    (In.m_bndCondExpansions),
            m_bndConditions        (In.m_bndConditions)
        {
            if(DeclareLinesSetCoeffPhys)
            {
                bool False = false;
                DisContField1DSharedPtr zero_line = boost::dynamic_pointer_cast<DisContField1D> (In.m_lines[0]);
                
                for(int n = 0; n < m_lines.num_elements(); ++n)
                {
                    m_lines[n] = MemoryManager<DisContField1D>::AllocateSharedPtr(*zero_line);
                }
                
                SetCoeffPhys();
            }
        }

        DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(LibUtilities::SessionReaderSharedPtr &pSession,
                                                                 const LibUtilities::BasisKey &HomoBasis_y,
																 const LibUtilities::BasisKey &HomoBasis_z,
																 const NekDouble lhom_y,
																 const NekDouble lhom_z,
																 bool useFFT,
																 SpatialDomains::MeshGraphSharedPtr &graph1D,
																 SpatialDomains::BoundaryConditions &bcs,
																 const int bc_loc):
            ExpList3DHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            int i,j,n,nel;
            bool True  = true; 
            bool False = false; 
            DisContField1DSharedPtr line_zero;

            //  
            m_lines[0] = line_zero = MemoryManager<DisContField1D>::AllocateSharedPtr(pSession,graph1D,bcs,bc_loc);

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_lines[0]->GetExpSize();

            for(j = 0; j < nel; ++j)
            {
                (*m_exp).push_back(m_lines[0]->GetExp(j));
            }
			
			int nylines = m_homogeneousBasis_y->GetNumPoints();
			int nzlines = m_homogeneousBasis_z->GetNumPoints();

            for(n = 1; n < nylines*nzlines; ++n)
            {
                m_lines[n] = MemoryManager<DisContField1D>::AllocateSharedPtr(pSession,graph1D,bcs,bc_loc);
                for(j = 0; j < nel; ++j)
                {
                    (*m_exp).push_back((*m_exp)[j]);
                }
            }            

            // Setup Default optimisation information. 
            nel = GetExpSize();

            m_globalOptParam = MemoryManager<NekOptimize::GlobalOptParam>
                ::AllocateSharedPtr(nel);
            
            SetCoeffPhys();

            SetupBoundaryConditions(HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,bcs);
        }

        DisContField3DHomogeneous2D::~DisContField3DHomogeneous2D()
        {
        }


        void DisContField3DHomogeneous2D::SetupBoundaryConditions(const LibUtilities::BasisKey &HomoBasis_y,
																  const LibUtilities::BasisKey &HomoBasis_z,
																  const NekDouble lhom_y,
																  const NekDouble lhom_z,
																  SpatialDomains::BoundaryConditions &bcs)
        {
            int i,j,n;
            
			// Setup an ExpList1DHomogeneous2D expansion for boundary
            // conditions and link to class declared in m_lines.
			
			int nlines = m_lines.num_elements();
			
			SpatialDomains::BoundaryRegionCollection  &bregions = bcs.GetBoundaryRegions();
			
			int nbnd = bregions.size();
			
			
			m_bndCondExpansions  = Array<OneD,MultiRegions::ExpListSharedPtr>(nbnd);
			
			Array<OneD, MultiRegions::ExpListSharedPtr> LinesBndCondExp(nlines);
			
			m_bndConditions = m_lines[0]->UpdateBndConditions();
            
            for(i = 0; i < nbnd; ++i)
            {
                for(n = 0; n < nlines; ++n)
                {
                    LinesBndCondExp[n] = m_lines[n]->UpdateBndCondExpansion(i);
                }
                
                m_bndCondExpansions[i] = MemoryManager<ExpList1DHomogeneous2D>::AllocateSharedPtr(m_session,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,m_useFFT,LinesBndCondExp);
                
            }
            
            EvaluateBoundaryConditions();
        }

        void DisContField3DHomogeneous2D::EvaluateBoundaryConditions(const NekDouble time)
        {
            int n,m;
			
			const Array<OneD, const NekDouble> y = m_homogeneousBasis_y->GetZ();
            const Array<OneD, const NekDouble> z = m_homogeneousBasis_z->GetZ();
			

            for(n = 0; n < m_nz; ++n)
            {
				for(m = 0; m < m_ny; ++m)
				{
					m_lines[m+(n*m_ny)]->EvaluateBoundaryConditions(time,0.5*m_lhom_y*(1.0+y[m]),0.5*m_lhom_z*(1.0+z[n]));
				}
            }
            
            // Fourier transform coefficient space boundary values
            for(n = 0; n < m_bndCondExpansions.num_elements(); ++n)
            {
                m_bndCondExpansions[n]->HomogeneousFwdTrans(m_bndCondExpansions[n]->GetCoeffs(),m_bndCondExpansions[n]->UpdateCoeffs());
            }    
        }
        
        void DisContField3DHomogeneous2D::v_HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray,
                    NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff)
        {
            v_HelmSolveDG(inarray, outarray, lambda, varLambda, varCoeff, 1);
        }

        void DisContField3DHomogeneous2D::v_HelmSolveDG(
                    const Array<OneD, const NekDouble> &inarray,
                    Array<OneD,       NekDouble> &outarray,
                    NekDouble lambda,
                    const Array<OneD, const NekDouble> &varLambda,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff,
                    NekDouble tau)
        {
            int n,m;
            int cnt = 0;
            int cnt1 = 0;
            int nhom_modes_y = m_homogeneousBasis_y->GetNumModes();
			int nhom_modes_z = m_homogeneousBasis_z->GetNumModes();
            NekDouble beta_y;
			NekDouble beta_z;
			NekDouble beta;
			
            Array<OneD, NekDouble> e_out;
            Array<OneD, NekDouble> fce(inarray.num_elements());

            // Fourier transform forcing function
			if(m_FourierSpace != eCoef)
			{
				HomogeneousFwdTrans(inarray,fce);
			}

            for(n = 0; n < nhom_modes_z; ++n)
            {
				for(m = 0; m < nhom_modes_y; ++m)
				{
					beta_z = 2*M_PI*(n/2)/m_lhom_z;
					beta_y = 2*M_PI*(m/2)/m_lhom_y;
					beta = beta_y*beta_y + beta_z*beta_z;
					
					m_lines[n]->HelmSolve(fce + cnt,
                                         e_out = outarray + cnt1,
                                         lambda + beta,varLambda,
                                         varCoeff);
                
					cnt  += m_lines[n]->GetTotPoints();
					cnt1 += m_lines[n]->GetNcoeffs();
				}
			}
        }
		
		void DisContField3DHomogeneous2D::v_EvaluateBoundaryConditions(const NekDouble time,const NekDouble x2_in, const NekDouble x3_in)
		{
			EvaluateBoundaryConditions(time);
		}
		
		const Array<OneD,const boost::shared_ptr<ExpList> > &DisContField3DHomogeneous2D::v_GetBndCondExpansions(void)
		{
			return GetBndCondExpansions();
		}
		
		const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> &DisContField3DHomogeneous2D::v_GetBndConditions()
		{
			return GetBndConditions();
		}
		
		boost::shared_ptr<ExpList> &DisContField3DHomogeneous2D::v_UpdateBndCondExpansion(int i)
		{
			return UpdateBndCondExpansion(i);
		}
		
		Array<OneD, SpatialDomains::BoundaryConditionShPtr> &DisContField3DHomogeneous2D::v_UpdateBndConditions()
		{
			return UpdateBndConditions();
		}


    } // end of namespace
} //end of namespace
