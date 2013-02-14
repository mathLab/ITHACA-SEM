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
#include <MultiRegions/DisContField1D.h>

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

        DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                                 const LibUtilities::BasisKey &HomoBasis_y,
																 const LibUtilities::BasisKey &HomoBasis_z,
                                                                 const NekDouble lhom_y,
																 const NekDouble lhom_z,
																 const bool useFFT,
																 const bool dealiasing):
            ExpList3DHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT,dealiasing),
            m_bndCondExpansions(),
            m_bndConditions()
        {
        }

        DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(const DisContField3DHomogeneous2D &In, const bool DeclareLinesSetCoeffPhys):
            ExpList3DHomogeneous2D (In,false),
            m_bndCondExpansions    (In.m_bndCondExpansions),
            m_bndConditions        (In.m_bndConditions)
        {
            if(DeclareLinesSetCoeffPhys)
            {
                DisContField1DSharedPtr zero_line = boost::dynamic_pointer_cast<DisContField1D> (In.m_lines[0]);
                
                for(int n = 0; n < m_lines.num_elements(); ++n)
                {
                    m_lines[n] = MemoryManager<DisContField1D>::AllocateSharedPtr(*zero_line);
                }
                
                SetCoeffPhys();
            }
        }

        DisContField3DHomogeneous2D::DisContField3DHomogeneous2D(const LibUtilities::SessionReaderSharedPtr &pSession,
                                                                 const LibUtilities::BasisKey &HomoBasis_y,
																 const LibUtilities::BasisKey &HomoBasis_z,
																 const NekDouble lhom_y,
																 const NekDouble lhom_z,
																 const bool useFFT,
																 const bool dealiasing,
																 const SpatialDomains::MeshGraphSharedPtr &graph1D,
																 const std::string &variable):
            ExpList3DHomogeneous2D(pSession,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,useFFT,dealiasing),
            m_bndCondExpansions(),
            m_bndConditions()
        {
            int i,n,nel;
            DisContField1DSharedPtr line_zero;
            SpatialDomains::BoundaryConditions bcs(pSession, graph1D);

            //  
            m_lines[0] = line_zero = MemoryManager<DisContField1D>::AllocateSharedPtr(pSession,graph1D,variable);

            m_exp = MemoryManager<StdRegions::StdExpansionVector>::AllocateSharedPtr();
            nel = m_lines[0]->GetExpSize();

            for(i = 0; i < nel; ++i)
            {
                (*m_exp).push_back(m_lines[0]->GetExp(i));
            }
			
			int nylines = m_homogeneousBasis_y->GetNumPoints();
			int nzlines = m_homogeneousBasis_z->GetNumPoints();

            for(n = 1; n < nylines*nzlines; ++n)
            {
                m_lines[n] = MemoryManager<DisContField1D>::AllocateSharedPtr(pSession,graph1D,variable);
                for(i = 0; i < nel; ++i)
                {
                    (*m_exp).push_back((*m_exp)[i]);
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
            int i,n;
            
			// Setup an ExpList1DHomogeneous2D expansion for boundary
            // conditions and link to class declared in m_lines.
			
			int nlines = m_lines.num_elements();
			
			const SpatialDomains::BoundaryRegionCollection  &bregions = bcs.GetBoundaryRegions();
			
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
                
                m_bndCondExpansions[i] = MemoryManager<ExpList1DHomogeneous2D>::AllocateSharedPtr(m_session,HomoBasis_y,HomoBasis_z,lhom_y,lhom_z,m_useFFT,false,LinesBndCondExp);
                
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
                const FlagList &flags,
                const StdRegions::ConstFactorMap &factors,
                const StdRegions::VarCoeffMap &varcoeff,
                const Array<OneD, const NekDouble> &dirForcing)
        {
            int n,m;
            int cnt = 0;
            int cnt1 = 0;
            int nhom_modes_y = m_homogeneousBasis_y->GetNumModes();
			int nhom_modes_z = m_homogeneousBasis_z->GetNumModes();
            NekDouble beta_y;
			NekDouble beta_z;
			StdRegions::ConstFactorMap new_factors;
			
            Array<OneD, NekDouble> e_out;
            Array<OneD, NekDouble> fce(inarray.num_elements());

            // Fourier transform forcing function
			if(m_WaveSpace)
			{
				fce = inarray;
			}
			else 
			{
				HomogeneousFwdTrans(inarray,fce);
			}

            for(n = 0; n < nhom_modes_z; ++n)
            {
				for(m = 0; m < nhom_modes_y; ++m)
				{
					beta_z = 2*M_PI*(n/2)/m_lhom_z;
					beta_y = 2*M_PI*(m/2)/m_lhom_y;
					new_factors = factors;
					new_factors[StdRegions::eFactorLambda] += beta_y*beta_y + beta_z*beta_z;
					
					m_lines[n]->HelmSolve(fce + cnt,
                                         e_out = outarray + cnt1,
                                         flags, new_factors, varcoeff, dirForcing);
                
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
		
		void DisContField3DHomogeneous2D::GetBoundaryToElmtMap(Array<OneD, int> &ElmtID, Array<OneD,int> &EdgeID)
        {
            
            if(m_BCtoElmMap.num_elements() == 0)
            {
                Array<OneD, int> ElmtID_tmp;
                Array<OneD, int> EdgeID_tmp;
				
                m_lines[0]->GetBoundaryToElmtMap(ElmtID_tmp,EdgeID_tmp);
                int nel_per_lines = m_lines[0]->GetExpSize();
                int nlines = m_lines.num_elements();
				
                int MapSize = ElmtID_tmp.num_elements();
				
                ElmtID = Array<OneD, int>(nlines*MapSize);
                EdgeID = Array<OneD, int>(nlines*MapSize);
                for(int i = 0; i < nlines; i++)
                {
                    for(int j = 0; j < MapSize; j++)
                    {
                        ElmtID[j+i*MapSize] = ElmtID_tmp[j] + i*nel_per_lines;
                        EdgeID[j+i*MapSize] = EdgeID_tmp[j];
                    }
                }
                m_BCtoElmMap = Array<OneD, int>(nlines*MapSize);
                m_BCtoEdgMap = Array<OneD, int>(nlines*MapSize);
				
                Vmath::Vcopy(nlines*MapSize,ElmtID,1,m_BCtoElmMap,1);
                Vmath::Vcopy(nlines*MapSize,EdgeID,1,m_BCtoEdgMap,1);
            }
            else 
            {
                int MapSize = m_BCtoElmMap.num_elements();
				
                ElmtID = Array<OneD, int>(MapSize);
                EdgeID = Array<OneD, int>(MapSize);
				
                Vmath::Vcopy(MapSize,m_BCtoElmMap,1,ElmtID,1);
                Vmath::Vcopy(MapSize,m_BCtoEdgMap,1,EdgeID,1);
            }			
        }


    } // end of namespace
} //end of namespace
